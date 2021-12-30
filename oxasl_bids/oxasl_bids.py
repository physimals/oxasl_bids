"""
OXASL_BIDS: Maps BIDS data sets on to oxasl/oxford_asl options
"""
import os.path as op
import os
import copy
import logging
import shutil

import numpy as np
import bids

from .mappings import get_oxasl_config_from_metadata

LOG = logging.getLogger(__name__)

# Catch-all exception type for BIDS data that cannot be mapped onto a 
# corresponding oxford_asl argument 
class IncompatabilityError(Exception):
    pass 

def get_command_line(options, prog="oxasl", extra_args=[]):
    """
    :param options: Dictionary of options derived from BIDS

    :return: Command string to run oxasl or oxford_asl
    """
    txt = prog + ' '
    for key,val in options.items():
        if key == "asl":
            key = "-i"
        elif key == "calib":
            key = "-c"
        elif key == "struct":
            key = "-s"
        elif key == "tes":
            continue # FIXME multi-TE support
        else:
            key = "--" + key

        if isinstance(val, list):
            val = ",".join([str(v) for v in val])

        if isinstance(val, bool):
            if val: 
                txt += f"{key} "
        elif val is not None:
            txt += f"{key} {val} "

    txt += " ".join(extra_args)
    return txt 

def _guarded_merge(common, specific):
    """
    Merge a set of oxasl options with a user-specified dictionary of common
    options. Common options override BIDS options (with a warning)
    
    :param common: User supplied dictionary
    :param specific: BIDS-derived options
    """
    merged = copy.deepcopy(specific)
    for aidx, arg_string in enumerate(common):
        arg_string = arg_string.strip('-')
        if arg_string.count('='):
            key = arg_string[:arg_string.index('=')]
            val = arg_string[arg_string.index('=')+1:]

            if key in specific:
                LOG.warn(f"Overwriting --{key} with value set in --common_args " + 
                f"({val} replaces {specific[key]}).")
            merged[key] = val 

        else: 
            merged[arg_string] = True

    return merged

def _get_asl_config(asl_file):
    """
    Identify the location of the ASL and calibration data
    and the ASL data order
    """
    options = {"asl" : op.abspath(asl_file.path)}
    options.update(get_oxasl_config_from_metadata(asl_file.get_metadata(), "asl"))
    metadata = asl_file.get_metadata()

    # Get the ASL context and interpret it
    ctx_filename = None
    for bids_file in asl_file.get_associations():
        if bids_file.entities["suffix"] == "aslcontext":
            ctx_filename = bids_file.path

    if ctx_filename is None:
        ctx_filename = asl_file.path[:asl_file.path.index(".nii")] + "context.tsv"
        ctx_filename = os.path.join(asl_file.dirname, ctx_filename)

    if not os.path.isfile(ctx_filename):
        raise IncompatabilityError("ASL context file not found")

    with open(ctx_filename) as f:
        ctx = [l.lower().strip() for l in f.readlines() if l.strip()]
    if ctx[0] != 'volume_type':
        raise IncompatabilityError("First line in ASL context is not volume_type")
    ctx = ctx[1:]
    if 'cbf' in ctx:
        raise IncompatabilityError("ASL context contains CBF images")
    if 'deltam' in ctx and ('label' in ctx or 'control' in ctx):
        raise IncompatabilityError("ASL sequence is mixed deltam and control/label")

    keys = {'control' : 0,
            'control(perev)': 0,
            'label': 1,
            'deltam': 2,
            'm0scan': 3}

    frame_codes = np.array([ keys[frame] for frame in ctx ])
    asl_frames = np.array([ i for i,code in enumerate(frame_codes) if code < 3 ])
    calib_frames = np.array([ i for i,code in enumerate(frame_codes) if code == 3 ])

    if frame_codes[asl_frames[0]] == 2:
        options["iaf"] = 'diff'
    elif frame_codes[asl_frames[0]] > frame_codes[asl_frames[1]]:
        options["iaf"] = 'tc'
    elif frame_codes[asl_frames[0]] < frame_codes[asl_frames[1]]:
        options["iaf"] = 'ct'
    else: 
        raise RuntimeError("Could not interpret label-control order")

    if calib_frames.size:
        # We have m0scan volumes in the ASL context - this suggests M0 is included in
        # ASL data, so check this and determine volume index
        if metadata.get('M0Type', 'included').lower() != "included":
            raise IncompatabilityError("JSON M0Type field not set to 'included', but m0scan in aslcontext")
        LOG.debug(f"Extracting M0 from {asl_file.filename}")
        options["calib"] = op.abspath(asl_file.path)
        options["calib_volumes"] = calib_frames
        options["asl_volumes"] = asl_frames
        options.update(get_oxasl_config_from_metadata(asl_file.get_metadata(), "calib"))
    else: 
        # No sign of m0scan volumes in ASL context - check M0 type is separate
        # and look for it in associated files
        if metadata.get('M0Type', 'separate').lower() != "separate":
            raise IncompatabilityError("JSON M0Type field not set to 'separate', but m0scan not in aslcontext")
        for bids_file in asl_file.get_associations():
            if bids_file.entities["suffix"] == "m0scan":
                options["calib"] = op.abspath(bids_file.path)
                options.update(get_oxasl_config_from_metadata(bids_file.get_metadata(), "calib"))
                LOG.debug(f"Found M0 in separate file: {bids_file.filename}")
    return options

def _get_struct_config(struct_file):
    return {"struct" : op.abspath(struct_file.path)}

def _get_calib_config(m0_file):
    # FIXME
    return {}

def _get_cblip_config(m0_file):
    ret = {}
    if m0_file.entities.get("datatype", None) == "fmap":
        ret["cblip"] = op.abspath(m0_file.path)
        ret.update(get_oxasl_config_from_metadata(m0_file.get_metadata(), "cblip", m0_file.get_image().shape))

    return ret

def _get_oxasl_config(asl_file, sess_files):
    """
    Build configuration options from JSON metadata
    and ASL context file
    """
    options = _get_asl_config(asl_file)
    for bids_file in asl_file.get_associations():
        if bids_file.entities["suffix"] == "m0scan":
            # Process associated M0 file (calib or cblip)
            pass

    for struct_file in sess_files["T1w"]:
        options.update(_get_struct_config(struct_file))

    for calib_file in sess_files["m0scan"]:
        if "calib" not in options:
            options.update(_get_calib_config(calib_file))
        options.update(_get_cblip_config(calib_file))

    return options

def get_bids_asl_files(bids_root):
    """
    :param bids_root: Path to root of BIDS dataset

    :return dict mapping subjects to sessions. Each session is a dict mapping
            file suffixes asl, m0scan, and T1w to lists of corresponding BIDSImageFile
            instances
    """
    dataset = bids.BIDSLayout(bids_root)
    data_files = {}
    for subj in dataset.get_subjects():
        data_files[subj] = {}
        sessions = dataset.get_sessions(subject=subj)
        if not sessions:
            sessions = [None]
        for sess in sessions:
            LOG.info("Getting ASL data for subject %s in session %s" % (subj, sess))
            data_files[subj][sess] = {}
            for suffix in ("asl", "m0scan", "T1w"):
                data_files[subj][sess][suffix] = []
                for fobj in dataset.get(subject=subj, session=sess, suffix=suffix):
                    if isinstance(fobj, bids.layout.models.BIDSImageFile):
                        LOG.debug("Found %s image: %s %s" % (suffix.upper(), fobj.filename, fobj.entities["suffix"]))
                        data_files[subj][sess][suffix].append(fobj)

    return data_files

def oxasl_config_from_bids(bids_root, common_options=None):
    """
    :param bids_root: Path to root of BIDS dataset
    :param common_options: Optional dictionary of oxasl options to add to BIDS derived options

    :return Sequence of OXASL configuration options, one for each ASL file found
            in the BIDS dataset
    """
    configs = []
    asl_files = get_bids_asl_files(bids_root)
    for subj, sessions in asl_files.items():
        for sess, sess_files in sessions.items():
            for asl_file in sess_files["asl"]:
                bids_options = _get_oxasl_config(asl_file, sess_files)
                if common_options:
                    bids_options = _guarded_merge(common_options, bids_options)
                LOG.debug("OXASL config for file %s" % asl_file.filename)
                LOG.debug(get_command_line(bids_options))
                configs.append(bids_options)
    return configs

def get_fslanat_command(config):
    struct_data = config.pop("struct")
    config["fslanat"] = os.path.basename(struct_data[:struct_data.index(".nii")]) + ".anat"
    return "fsl_anat -i %s\n" % struct_data

OXASL_OUTPUT_MAPPING = [
    # FIXME .nii.gz extension?
    ("native_space/perfusion_calib.nii.gz", "_CBFmap.nii.gz"),
    ("native_space/aCBV_calib.nii.gz", "aCBVmap.nii.gz"),
    ("native_space/arrival.nii.gz", "ATTmap.nii.gz"),
    ("native_space/perfusion_var_calib.nii.gz", "CBFvar.nii.gz"),
    ("native_space/aCBV_var_calib.nii.gz", "aCBVvar.nii.gz"),
    ("native_space/arrival_var.nii.gz", "ATTvar.nii.gz"),
    ("native_space/mask.nii.gz", "mask.nii.gz"),
    # FIXME normalised / relative CBF, sensitivity map, WM/GM masks, 
    # cortical/cerebral maps
    # Transformation matrices
    # GM/WM mean values
]

def copy_bids_dataset(bidsdir, output_dir, subject=None, session=None):
    raise NotImplementedError()

def bids_filename(suffix, subject, session):
    if session:
        return "sub-%s_ses-%s_%s" % (subject, session, suffix)
    else:
        return "sub-%s_%s" % (subject, suffix)

def oxford_asl_to_bids(oxford_asl_dir, bidsdir, subject, session=None, bids_output_dir=None):
    """
    Convert oxford_asl output to BIDS format
    
    :param oxasl_dir: Oxford_asl output folder
    :param output_dir: Destination BIDS output
    :param bids_src: Source BIDS dataset. If not specified we will assume the output
                     dir already contains a BIDS dataset which we are merging into
    
    """
    if bids_output_dir:
        # We are creating a separate output dataset, so we need to start by copying
        # the BIDS source dataset
        copy_bids_dataset(bidsdir, bids_output_dir, subject, session)
    else:
        bids_output_dir = bidsdir

    base_dir = os.path.join(bids_output_dir, "derivatives", "oxford_asl", "sub-%s" % subject)
    if session:
        base_dir = os.path.join(base_dir, "ses-%s" % session)
    os.makedirs(base_dir, exist_ok=True)

    for (src, dest) in OXASL_OUTPUT_MAPPING:
        dest =  os.path.join(base_dir, bids_filename(dest, subject, session))
        src = os.path.join(oxford_asl_dir, src)
        if op.exists(src):
            shutil.copy(src, dest)
        else:
            LOG.warn("Oxford_asl output file not found: %s" % src)
