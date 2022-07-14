"""
OXASL_BIDS: Miscellaneous utilities
"""
import json

import nibabel as nib

def bids_filename(suffix, subject, session, labeldict=None):
    """
    Get a BIDS style filename
    
    :param suffix: The image type suffix, e.g. T1w
    :param subject: Subject ID
    :param session: Session ID
    :param labeldict: Dictionary of labels to be embedded in the filename
    """
    fname = f"sub-{subject}"
    if session:
        fname += f"_ses-%{session}"
    if labeldict:
        for key, value in labeldict.items():
            fname += f"_{key}-{value}"
    fname += f"_{suffix}"
    return fname

def load_img(fname):
    """
    Load an image

    :return: Tuple of nii structure, JSON metadata dictionary
    """
    nii = nib.load(fname)
    json_filename = fname[:fname.index(".nii")] + ".json"
    with open(json_filename, "r") as f:
        metadata = json.load(f)
    metadata["img_shape"] = nii.shape

    return nii, metadata
