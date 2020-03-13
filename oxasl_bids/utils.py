import warnings

import nibabel
import numpy as np 

# Local function to read out an FSL-specific affine matrix from an image
def _getFSLspace(imgPth):
    obj = nibabel.load(imgPth)
    if obj.header['dim'][0] < 3:
        raise RuntimeError("Volume has less than 3 dimensions" + \
                "cannot resolve space")

    sform = obj.affine
    det = np.linalg.det(sform[0:4, 0:4])
    ret = np.identity(4)
    pixdim = obj.header['pixdim'][1:4]
    for d in range(3):
        ret[d,d] = pixdim[d]

    # Check the xyzt field to find the spatial units. 
    xyzt =str(obj.header['xyzt_units'])
    if xyzt == '01': 
        multi = 1000
    elif xyzt == '10':
        multi = 1 
    elif xyzt =='11':
        multi = 1e-3
    else: 
        multi = 1
        warnings.warn("Assuming mm units for transform")

    if det > 0:
        ret[0,0] = -pixdim[0]
        ret[0,3] = (obj.header['dim'][1] - 1) * pixdim[0]

    ret = ret * multi
    ret[3,3] = 1
    return ret


def _FLIRT_to_world(source, reference, transform):
    """Adjust a FLIRT transformation matrix into a true world-world 
    transform. Required as FSL matrices are encoded in a specific form 
    such that they can only be applied alongside the requisite images (extra
    information is required from those images). With thanks to Martin Craig
    and Tim Coalson. See: https://github.com/Washington-University/workbench/blob/9c34187281066519e78841e29dc14bef504776df/src/Nifti/NiftiHeader.cxx#L168 
    https://github.com/Washington-University/workbench/blob/335ad0c910ca9812907ea92ba0e5563225eee1e6/src/Files/AffineFile.cxx#L144

    Args: 
        source: path to source image, the image to be deformed 
        reference: path to reference image, the target of the transform
        transform: affine matrix produced by FLIRT from src to ref 

    Returns: 
        complete transformation matrix between the two. 
    """

    # 'Space' refers to conversion of world -> scaled voxel coordinates, 
    # with no shift for origin 
    Ss = _getFSLspace(source)
    Rs = _getFSLspace(reference)

    # 'Affine' refers to the voxel -> world conversion matrix 
    Ra = nibabel.load(reference).affine
    Sa = nibabel.load(source).affine

    # Work backwards (R->L) to interpret these 
    out = Ra @ np.linalg.inv(Rs) @ transform @ Ss @ np.linalg.inv(Sa)
    return out


def _world_to_FLIRT(source, reference, transform):
    """
    Inverse of _FLIRT_to_world: convert a world-world transform to FLIRT matrix
    """

    # 'Space' refers to conversion of world -> scaled voxel coordinates, 
    # with no shift for origin 
    Ss = _getFSLspace(source)
    Rs = _getFSLspace(reference)

    # 'Affine' refers to the voxel -> world conversion matrix 
    Ra = nibabel.load(reference).affine
    Sa = nibabel.load(source).affine

    # Work backwards (R->L) to interpret these 
    out = Rs @ np.linalg.inv(Ra) @ transform @ Sa @ np.linalg.inv(Ss)
    return out

class IncompatabilityError(Exception):
    pass 