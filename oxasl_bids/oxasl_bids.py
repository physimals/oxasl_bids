"""
OXASL_BIDS: Maps BIDS data sets on to oxasl options and oxasl output back to BIDS
"""
import os.path as op
import os
import copy
import logging
import shutil

import bids

from . import utils
from .mappings import oxasl_config_from_metadata

LOG = logging.getLogger(__name__)

OXASL_OUTPUT_MAPPING = {
    "perfusion_calib" : "CBFmap",
    "aCBV_calib" : "aCBVmap",
    "arrival" : "ATTmap",
    "perfusion_var_calib" : "CBFvar",
    "aCBV_var_calib" : "aCBVvar",
    "arrival_var" : "ATTvar",
    "mask" : "mask"
    # FIXME normalised / relative CBF, sensitivity map, WM/GM masks, 
    # cortical/cerebral maps
    # Transformation matrices
    # GM/WM mean values
    # Structural / std space outputs
    # PVC outputs
}

def oxasl_config_from_bids(bids_root, common_options=None):
    """
    Get OXASL configuration options from a BIDS data set

    :param bids_root: Path to root of BIDS dataset
    :param common_options: Optional dictionary of oxasl options to add to BIDS derived options

    :return Sequence of OXASL configuration options, one for each ASL file found
            in the BIDS dataset
    """
    configs = []
    bids_sessions = _get_bids_sessions(bids_root)
    for subjid, sessions in bids_sessions.items():
        for sessid, sess_files in sessions.items():
            for asl_file in sess_files["asl"]:
                bids_options = _get_oxasl_config(asl_file, sess_files)
                if common_options:
                    bids_options = _guarded_merge(common_options, bids_options)
                if "output" not in bids_options:
                    output_dir = f"sub-{subjid}"
                    if sessid:
                        output_dir += f"sess-{sessid}"
                    bids_options["output"] = output_dir

                LOG.debug("OXASL config for file %s" % asl_file.filename)
                LOG.debug(get_oxasl_command_line(bids_options))
                configs.append({"options" : bids_options, "subject" : subjid, "session" : sessid})
    return configs

def oxasl_output_to_bids(oxasl_dir, bidsdir, subject, session=None, bids_output_dir=None):
    """
    Convert oxasl output to BIDS format
    
    :param oxasl_dir: Oxasl output folder
    :param bidsdir: Source BIDS dataset. If not specified we will assume the output
                    dir already contains a BIDS dataset which we are merging into
    :param subject: Subject ID for this oxasl run
    :param session: Session ID for this oxasl run. If not specified assume only
                    one session present in dataset
    :param output_dir: Destination BIDS output. If not specified, merge with source BIDS
                       data as a derivative
    """
    if bids_output_dir:
        # We are creating a separate output dataset, so we need to start by copying
        # the BIDS source dataset
        # FIXME this isn't quite right we should create a separate data set with optional
        # sourcedata and maybe it could already exist and include multiple subjects...
        #_copy_bids_dataset(bidsdir, bids_output_dir, subject, session)
        bids_output_subdir = "derivatives"
    else:
        bids_output_dir = bidsdir
        bids_output_subdir = "derivatives"

    # FIXME, check that subject exists in bidsdir and session too (or if no session that
    # source dataset only has single session for this subject)
    base_dir = os.path.join(bids_output_dir, bids_output_subdir, "oxasl", "sub-%s" % subject)
    if session:
        base_dir = os.path.join(base_dir, "ses-%s" % session)
    os.makedirs(base_dir, exist_ok=True)

    for space in ("native", "std", "struct"):
        srcdir = os.path.join(oxasl_dir, f"{space}_space")
        LOG.info(f"Looking for {space} space output")
        if not os.path.isdir(srcdir): continue
        for src, dest in OXASL_OUTPUT_MAPPING.items():
            srcpath, ext = _getnii(os.path.join(srcdir, src))
            if srcpath:
                destpath =  os.path.join(base_dir, utils.bids_filename(dest + ext, subject, session, {"space" : space}))
                LOG.info(f"Copying {srcpath} to {destpath}")
                shutil.copy(srcpath, destpath)
            else:
                LOG.warn("Oxasl output file not found: %s" % src)

def get_output_as_bids_command(args, config):
    """
    Get the command to convert the output of oxasl to a BIDS data set

    :param args: Command line arguments
    :param config: oxasl session config
    """
    options = config["options"]
    cmdline = f"oxasl_bids bidsout --bidsdir {os.path.abspath(args.bidsdir)} --oxasl-output {options['output']} --subject {config['subject']}"
    if config['session']:
        cmdline += f" --session {config['session']}"
    if args.bids_output:
        cmdline += f" --bids-output {os.path.abspath(args.bids_output)}"
    if args.merge_source:
        cmdline += f" --merge-source"
    return cmdline

def get_fslanat_command(options):
    """
    Update OXASL configuration to use a separately run FSL_ANAT command on the structural
    image rather than passing it to OXASL directly.

    :param options: OXASL configuration dict. Structural image option will be removed
                   and FSL_ANAT option added
    :return: Command line to run FSL_ANAT on the structural data from the OXASL configuration
    """
    struct_data = options.pop("struct")
    struct_name = os.path.basename(struct_data[:struct_data.index(".nii")])
    options["fslanat"] = struct_name + ".anat"
    return f"fsl_anat -i {struct_data} -o {options['fslanat']}\n"

def get_oxasl_command_line(options, extra_args=[]):
    """
    :param options: Dictionary of options derived from BIDS

    :return: Command string to run oxasl
    """
    txt = 'oxasl '
    for key,val in options.items():
        if key == "asl":
            key = "-i"
        elif key == "calib":
            key = "-c"
        elif key == "struct":
            key = "-s"
        elif key == "output":
            key = "-o"
        elif key == "tes":
            continue # FIXME multi-TE support
        else:
            key = "--" + key

        if isinstance(val, list):
            if all([isinstance(v, int) for v in val]):
                val = ",".join(["%i" % v for v in val])
            else:
                val = ",".join(["%.3g" % v for v in val])
        elif isinstance(val, float):
            val = "%.4g" % val

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
                LOG.warn(f"Overwriting --{key} with value set in --common_args ({val} replaces {specific[key]})")
            merged[key] = val 

        else: 
            merged[arg_string] = True

    return merged

def _parse_aslcontext(ctx_filename):
    with open(ctx_filename) as f:
        ctx = [l.lower().strip() for l in f.readlines() if l.strip()]
    if ctx[0] != 'volume_type':
        raise utils.IncompatabilityError("First line in ASL context is not volume_type")
    ctx = ctx[1:]
    if 'cbf' in ctx:
        raise utils.IncompatabilityError("ASL context contains CBF images")
    if 'deltam' in ctx and ('label' in ctx or 'control' in ctx):
        raise utils.IncompatabilityError("ASL sequence is mixed deltam and control/label")

    keys = {'control' : 0,
            'control(perev)': 0,
            'label': 1,
            'deltam': 2,
            'm0scan': 3}

    frame_codes = [ keys[frame] for frame in ctx ]
    asl_frames = [ i for i,code in enumerate(frame_codes) if code < 3 ]
    calib_frames = [ i for i,code in enumerate(frame_codes) if code == 3 ]
    asl_frame_codes = [frame_codes[i] for i in asl_frames]

    if not asl_frames:
        raise utils.IncompatabilityError("No ASL data (label/control or deltam found in ASL data file")
    elif len(asl_frames) == 1 and frame_codes[asl_frames[0]] != 2:
        raise utils.IncompatabilityError("Only one ASL volume found in ASL data file and it was not a deltam image")

    return asl_frames, asl_frame_codes, calib_frames

def _get_iaf(asl_frame_codes):
    if asl_frame_codes[0] == 2:
        return 'diff'
    else:
        label_idxs = [idx for idx, frame_type in enumerate(asl_frame_codes) if frame_type == 1]
        control_idxs = [idx for idx, frame_type in enumerate(asl_frame_codes) if frame_type == 2]
        if all([idx % 2 == 0 for idx in label_idxs]) and all([idx % 2 == 1 for idx in control_idxs]):
            return 'tc'
        elif all([idx % 2 == 1 for idx in label_idxs]) and all([idx % 2 == 0 for idx in control_idxs]):
            return 'ct'
        else:
            # FIXME detect other kinds of ordering, e.g. TTTTTCCCCC etc
            raise utils.IncompatabilityError("Could not interpret label-control order in aslcontext.tsv")

def _get_asl_config(asl_file):
    """
    Extract relevant oxasl configuration from a BIDSImageFile containing ASL data
    """
    options = {"asl" : op.abspath(asl_file.path)}
    options.update(oxasl_config_from_metadata(asl_file.get_metadata(), "asl"))
    metadata = asl_file.get_metadata()

    # Get the ASL context and interpret it. This is what tells us the ordering of
    # label/control image, or if the data is already differenced
    ctx_filename = None
    for bids_file in asl_file.get_associations():
        if bids_file.entities["suffix"] == "aslcontext":
            ctx_filename = bids_file.path

    if ctx_filename is None:
        ctx_filename = asl_file.path[:asl_file.path.index(".nii")] + "context.tsv"
        ctx_filename = os.path.join(asl_file.dirname, ctx_filename)

    if not os.path.isfile(ctx_filename):
        raise utils.IncompatabilityError("ASL context file not found")

    asl_frames, asl_frame_codes, calib_frames = _parse_aslcontext(ctx_filename)
    options["iaf"] = _get_iaf(asl_frame_codes)
    if options["iaf"] in ("ct", "tc"):
        # With TC or CT pairs the timings will be repeated. We don't care about the order since
        # no valid ASL sequence will have different timings for tag and control.
        ttype = 'plds' if options.get('casl', False) else "tis"
        options[ttype] = [options[ttype][idx] for idx in range(0, len(options[ttype]), 2)]

    if calib_frames:
        # We have m0scan volumes in the ASL context - this suggests M0 is included in
        # ASL data, so check this and determine volume index. Calibration image options
        # (e.g. TR, TE) are then to be derived from this image
        if metadata.get('M0Type', 'included').lower() != "included":
            raise utils.IncompatabilityError("JSON M0Type field not set to 'included', but m0scan in aslcontext")
        LOG.debug(f"Extracting M0 from {asl_file.filename}")
        options["calib"] = op.abspath(asl_file.path)
        options["calib_volumes"] = calib_frames
        options["asl_volumes"] = asl_frames
        options.update(oxasl_config_from_metadata(asl_file.get_metadata(), "calib"))
    else: 
        # No sign of m0scan volumes in ASL context - check M0 type is separate
        # and look for it in associated files
        if metadata.get('M0Type', 'separate').lower() != "separate":
            raise utils.IncompatabilityError("JSON M0Type field not set to 'separate', but m0scan not in aslcontext")
        for bids_file in asl_file.get_associations():
            if bids_file.entities["suffix"] == "m0scan":
                LOG.debug(f"Found M0 in separate file referenced from ASL data: {bids_file.filename}")
                options["calib"] = op.abspath(bids_file.path)
                options.update(oxasl_config_from_metadata(bids_file.get_metadata(), "calib"))
    return options

def _get_struct_config(struct_file):
    return {"struct" : op.abspath(struct_file.path)}

def _get_calib_config(m0_file):
    """
    Check if M0 file is a calibration image and if so add configuration
    to the options

    FIXME how to check if actual calib file properly
    """
    ret = {}
    if m0_file.entities.get("datatype", None) != "fmap":
        ret["calib"] = op.abspath(m0_file.path)
        ret.update(oxasl_config_from_metadata(m0_file.get_metadata(), "calib", img_shape=m0_file.get_image().shape))

def _get_cblip_config(m0_file):
    """
    Check if M0 file is a blipped calibration image and if so add configuration
    to the options
    """
    ret = {}
    if m0_file.entities.get("datatype", None) == "fmap":
        ret["cblip"] = op.abspath(m0_file.path)
        ret.update(oxasl_config_from_metadata(m0_file.get_metadata(), "cblip", img_shape=m0_file.get_image().shape))

    return ret

def _get_oxasl_config(asl_file, sess_files):
    """
    Build configuration options from JSON metadata and ASL context file

    :param asl_file: BIDSImageFile containing the ASL data
    :param sess_files: Dictionary mapping file type to BIDSImageFile for 
                       other relevant files in the session (e.g. T1, m0)

    :return: Dictionary of option name, option value
    """
    # Extract as much configuration as we can from the ASL data. This may
    # include calibration images
    options = _get_asl_config(asl_file)

    # Look for structural data
    if sess_files["T1w"]:
        options.update(_get_struct_config(sess_files["T1w"][0]))
        if len(sess_files["T1w"]) > 1:
            LOG.warn("Multiple T1w structural images found for ASL file %{asl_file} - using first")

    # There may be multiple M0 scan files, e.g. if we have a 'blipped' M0 image. Note that
    # we only use calibration files if we didn't get this info from the ASL data file
    for calib_file in sess_files["m0scan"]:
        if "calib" not in options:
            options.update(_get_calib_config(calib_file))
        if "cblip" not in options:
            options.update(_get_cblip_config(calib_file))

    # Need to fix PLDs a bit 
    return options

def _get_bids_sessions(bids_root):
    """
    Get structure describing all sessions found in a BIDS dataset

    :param bids_root: Path to root of BIDS dataset

    :return dict mapping subjects to a group of sessions. The session group is a dict mapping
            session ID to session dict. Each session is a dict mapping file suffixes asl, m0scan, 
            and T1w to lists of corresponding BIDSImageFile instances
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

def _copy_bids_dataset(bidsdir, output_dir, subject=None, session=None):
    raise NotImplementedError()

def _getnii(fname):
    for ext in (".nii.gz", ".nii"):
        LOG.info(f"Looking for {fname + ext}")
        if os.path.exists(fname + ext):
            return fname + ext, ext
    return None, None
