"""
Code relating to the structural image processing

Moved here temporarily while we decide what to do with it
"""
import argparse
import glob
import copy
import os
import os.path as op
import functools
import multiprocessing
import warnings

import nibabel
import numpy as np

import toblerone
from scipy.ndimage import affine_transform

def get_alignment_path(target, asl_dir):
    if target == 'anat':
        anat_dir = asl_dir.replace('asl', 'anat')
        anat_path = glob.glob(op.join(anat_dir, '*_T1w.nii.gz'))
        if len(anat_path) > 1: 
            warnings.warn("Found more than one T1 image")
        anat_path = anat_path[0]
        return anat_path
    else: 
        raise RuntimeError("Not yet implemented")

def prepare_alignment_target(target_path, asl_path):
    target = nibabel.load(target_path)
    asl = nibabel.load(asl_path)
    asl_vox = asl.header['pixdim'][1:4]
    target_vox = target.header['pixdim'][1:4]

    factor = asl_vox / target_vox
    dims = np.ceil(target.shape / factor).astype(np.int16) 
    new_affine = copy.deepcopy(target.affine)
    for idx in range(3):
        new_affine[:,idx] *= factor[idx]
    
    target_img = target.get_fdata()
    min_max = (target_img.min(), target_img.max())
    asl2target = np.linalg.inv(target.affine) @ new_affine
    vol = affine_transform(target_img, asl2target, 
        output_shape=dims, mode='constant')
    vol[vol < min_max[0]] = min_max[0]
    vol[vol > min_max[1]] = min_max[1]
    align_img = nibabel.Nifti2Image(vol.astype(np.float32), affine=new_affine)

    return align_img
def run_fsl_anat(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--bidsdir', required=True, type=str)
    parser.add_argument('--cores', default=1, type=int)
    # parser.add_argument()
    args = dict(vars(parser.parse_args(argv)))
    cores = args['cores']
    bids_root = args['bidsdir']

    jobs = [] 
    for anat_dir in walk_modality_dirs(bids_root, 'anat'):
        anat = glob.glob(op.join(anat_dir, 'sub-*_T1w.nii*'))
        if isinstance(anat, list):
            anat = anat[0]

        stub = op.split(anat)[1]
        outdir = derivatives_dir(anat_dir)
        os.makedirs(outdir, exist_ok=True)
        outname = op.join(outdir, stub[:stub.index('.nii')] + '.anat')
        jobs.append({'struct': anat, 'out': outname})
        

    worker = functools.partial(starstarmap, toblerone.fsl_fs_anat)
    if cores > 1:
        with multiprocessing.Pool(cores) as p: 
            p.map(worker, jobs)
    else: 
        for job in jobs: 
            worker(job)

def starstarmap(func, arg_dict):
    return func(**arg_dict)
