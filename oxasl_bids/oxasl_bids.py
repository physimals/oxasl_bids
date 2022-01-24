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

from . import utils
from .mappings import get_oxasl_config_from_metadata

LOG = logging.getLogger(__name__)

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
                LOG.warn(f"Overwriting --{key} with value set in --common_args " + 
                f"({val} replaces {specific[key]}).")
            merged[key] = val 

        else: 
            merged[arg_string] = True

    return merged

def _get_asl_config(asl_file):
    """
    Extract relevant oxasl configuration from a BIDSImageFile containing ASL data
    """
    options = {"asl" : op.abspath(asl_file.path)}
    options.update(get_oxasl_config_from_metadata(asl_file.get_metadata(), "asl"))
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
    if not asl_frames:
        raise utils.IncompatabilityError("No ASL data (label/control or deltam found in ASL data file")
    elif len(asl_frames) == 1 and frame_codes[asl_frames[0]] != 2:
        raise utils.IncompatabilityError("Only one ASL volume found in ASL data file and it was not a deltam image")

    if frame_codes[asl_frames[0]] == 2:
        options["iaf"] = 'diff'
    else:
        label_idxs = [idx for idx, frame_type in enumerate(asl_frames) if frame_type == 1]
        control_idxs = [idx for idx, frame_type in enumerate(asl_frames) if frame_type == 2]
        if all([idx % 2 == 0 for idx in label_idxs]) and all([idx % 2 == 1 for idx in control_idxs]):
            options["iaf"] = 'tc'
        elif all([idx % 2 == 1 for idx in label_idxs]) and all([idx % 2 == 0 for idx in control_idxs]):
            options["iaf"] = 'ct'
        else:
            # FIXME detect other kinds of ordering, e.g. TTTTTCCCCC etc
            raise utils.IncompatabilityError("Could not interpret label-control order in aslcontext.tsv")

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
        options.update(get_oxasl_config_from_metadata(asl_file.get_metadata(), "calib"))
    else: 
        # No sign of m0scan volumes in ASL context - check M0 type is separate
        # and look for it in associated files
        if metadata.get('M0Type', 'separate').lower() != "separate":
            raise utils.IncompatabilityError("JSON M0Type field not set to 'separate', but m0scan not in aslcontext")
        for bids_file in asl_file.get_associations():
            if bids_file.entities["suffix"] == "m0scan":
                LOG.debug(f"Found M0 in separate file referenced from ASL data: {bids_file.filename}")
                options["calib"] = op.abspath(bids_file.path)
                options.update(get_oxasl_config_from_metadata(bids_file.get_metadata(), "calib"))
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
        ret.update(get_oxasl_config_from_metadata(m0_file.get_metadata(), "calib", m0_file.get_image().shape))

def _get_cblip_config(m0_file):
    """
    Check if M0 file is a blipped calibration image and if so add configuration
    to the options
    """
    ret = {}
    if m0_file.entities.get("datatype", None) == "fmap":
        ret["cblip"] = op.abspath(m0_file.path)
        ret.update(get_oxasl_config_from_metadata(m0_file.get_metadata(), "cblip", m0_file.get_image().shape))

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

def get_bids_sessions(bids_root):
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

def oxasl_config_from_bids(bids_root, common_options=None):
    """
    Get OXASL configuration options from a BIDS data set

    :param bids_root: Path to root of BIDS dataset
    :param common_options: Optional dictionary of oxasl options to add to BIDS derived options

    :return Sequence of OXASL configuration options, one for each ASL file found
            in the BIDS dataset
    """
    configs = []
    bids_sessions = get_bids_sessions(bids_root)
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
                LOG.debug(get_command_line(bids_options))
                configs.append({"options" : bids_options, "subject" : subjid, "session" : sessid})
    return configs

def get_output_as_bids_command(args, config):
    """
    Get the command to convert the output of oxford_asl to a BIDS data set

    :param args: Command line arguments
    :param config: oxford_asl session config
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
    return "fsl_anat -i %s -o %s\n" % (struct_data, options["fslanat"])

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

def copy_bids_dataset(bidsdir, output_dir, subject=None, session=None):
    raise NotImplementedError()

def bids_filename(suffix, subject, session, labeldict=None):
    fname = f"sub-{subject}"
    if session:
        fname += f"_ses-%{session}"
    if labeldict:
        for key, value in labeldict.items():
            fname += f"_{key}-{value}"
    fname += f"_{suffix}"
    return fname

def oxford_asl_to_bids(oxford_asl_dir, bidsdir, subject, session=None, bids_output_dir=None):
    """
    Convert oxford_asl output to BIDS format
    
    :param oxasl_dir: Oxford_asl output folder
    :param bidsdir: Source BIDS dataset. If not specified we will assume the output
                    dir already contains a BIDS dataset which we are merging into
    :param subject: Subject ID for this oxford_asl run
    :param session: Session ID for this oxford_asl run. If not specified assume only
                    one session present in dataset
    :param output_dir: Destination BIDS output. If not specified, merge with source BIDS
                       data as a derivative
    """
    if bids_output_dir:
        # We are creating a separate output dataset, so we need to start by copying
        # the BIDS source dataset
        # FIXME this isn't quite right we should create a separate data set with optional
        # sourcedata and maybe it could already exist and include multiple subjects...
        #copy_bids_dataset(bidsdir, bids_output_dir, subject, session)
        bids_output_subdir = "derivatives"
    else:
        bids_output_dir = bidsdir
        bids_output_subdir = "derivatives"

    # FIXME, check that subject exists in bidsdir and session too (or if no session that
    # source dataset only has single session for this subject)
    base_dir = os.path.join(bids_output_dir, bids_output_subdir, "oxford_asl", "sub-%s" % subject)
    if session:
        base_dir = os.path.join(base_dir, "ses-%s" % session)
    os.makedirs(base_dir, exist_ok=True)

    for space in ("native", "std", "struct"):
        srcdir = os.path.join(oxford_asl_dir, f"{space}_space")
        LOG.info(f"Looking for {space} space output")
        if not os.path.isdir(srcdir): continue
        for src, dest in OXASL_OUTPUT_MAPPING.items():
            srcpath, ext = _getnii(os.path.join(srcdir, src))
            if srcpath:
                destpath =  os.path.join(base_dir, bids_filename(dest + ext, subject, session, {"space" : space}))
                LOG.info(f"Copying {srcpath} to {destpath}")
                shutil.copy(srcpath, destpath)
            else:
                LOG.warn("Oxford_asl output file not found: %s" % src)

def _getnii(fname):
    for ext in (".nii.gz", ".nii"):
        LOG.info(f"Looking for {fname + ext}")
        if os.path.exists(fname + ext):
            return fname + ext, ext
    return None, None

# def wipe_dir(path):
#     shutil.rmtree(path)
#     os.makedirs(path, exist_ok=True)

# def prepare_config_files(argv): 
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--bidsdir', required=True, type=str)
#     parser.add_argument('--common_args', required=False, 
#         type=str, nargs=argparse.REMAINDER)
#     parser.add_argument('--align', type=str)
#     parser.add_argument('--overwrite', action='store_true')
#     parser.add_argument('--fsl_anat', action='store_true')

#     args = dict(vars(parser.parse_args(argv)))
#     bids_root = args.pop('bidsdir')
#     align_spc = args.pop('align')
#     fsl_anat = args.pop('fsl_anat')
#     overwrite = args.pop('overwrite')
#     common_args = args['common_args']

#     if align_spc and (align_spc != 'anat'):
#         raise RuntimeError("Not implemented yet")
#     if align_spc and (not fsl_anat):
#         raise RuntimeError("--align must be used with --fsl_anat") 

            # oxdir = oxasl_dir(asl_dir, asl_file)
            # if overwrite: 
            #     wipe_dir(oxdir)
            # elif os.listdir(oxdir):
            #     raise RuntimeError(f"Oxasl output directory {oxdir} is not empty. Use --overwrite option.")
            # config_dir = configuration_dir(asl_dir, asl_file)
            # outpath =  op.join(config_dir, 'oxasl_config.txt')
            # rel_path_root = op.abspath(op.join(oxdir, '..'))
            # outstring = _dump_to_string(bids_options, rel_path_root)
            # outstring += "--output={} ".format(oxdir)

            # if fsl_anat: 
            #     anatdir = asl_dir.replace('sourcedata', 'derivatives')
            #     anatdir = anatdir.replace('asl', 'anat')
            #     anatdirs = glob.glob(op.join(anatdir, 'sub-*T1w.anat'))
            #     if len(anatdirs) > 1:
            #         raise RuntimeError(f"Found multiple fsl_anat dirs in {anatdir}")
            #     elif not len(anatdirs):
            #         raise RuntimeError(f"Did not find fsl_anat dir in {anatdir}")
            #     fslanat = op.relpath(anatdirs[0], rel_path_root)
            #     outstring += "--fslanat={} ".format(fslanat)

            # if align_spc and fsl_anat: 
            #     t1 = op.join(anatdirs[0], 'T1_biascorr_brain.nii.gz')
            #     align_img = prepare_alignment_target(t1, op.join(asl_dir, asl))
            #     align_path = op.join(config_dir, 'oxasl_bids_common_space.nii.gz')
            #     nibabel.save(align_img, align_path)
            #     align_mat = utils._world_to_FLIRT(t1, align_path, np.eye(4))
            #     align_mat_path = op.join(config_dir, 'oxasl_bids_common_space_flirt.txt')
            #     rel_mat_path, rel_algn_path = [ op.relpath(p, rel_path_root) 
            #         for p in [ align_mat_path, align_path ] ] 
            #     np.savetxt(align_mat_path, align_mat)
            #     outstring += f' --output-custom={rel_algn_path} --output-custom-mat={rel_mat_path}'

            # file_count += 1     
            # outstring += ' --overwrite'
            # with open(outpath, 'w') as f: 
            #     f.write(outstring)

    #     except utils.IncompatabilityError as e:
    #         print(f"Warning: incompatible parameters found for {asl} - skipping")
    #         print(e)
    #         continue
    #     except FileNotFoundError as e: 
    #         print(f"Warning: missing file for {asl} - skipping")
    #         print(e)
    #         continue
    #     except Exception as e: 
    #         raise e 

    # print(f"Wrote {file_count} configuration files in {op.join(bids_root, 'derivatives')}")


# def __run_oxasl_worker(at_dir, cmd_path):
#     os.chdir(at_dir)
#     cmd = open(cmd_path, 'r').read()
#     return subprocess.run(cmd, shell=True)

# def run_config_files(argv):
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--bidsdir', required=True, type=str)
#     parser.add_argument('--cores', default=1, type=int)
#     args = dict(vars(parser.parse_args(argv)))
#     bids_root = args['bidsdir']
#     cores = args['cores']

#     jobs = []  
#     for derivs_dir in walk_modality_dirs(bids_root, 'perf', datatype_dir="derivatives"):
#         # Look for the oxasl_directory
#         all_dirs = glob.glob(op.join(derivs_dir, 'sub-*_oxasl'))
#         fltr = re.compile('sub-\d*.*_oxasl')
#         oxasl_dirs = [ p for p in all_dirs if fltr.match(op.split(p)[1]) ]
#         for oxdir in oxasl_dirs:
#             config = op.join(oxdir, 'config', 'oxasl_config.txt')
#             if not op.exists(config):
#                 LOG.warn(f"Expected to find a configuration file ({config}) in {oxdir}.")
#                 continue
#             else: 
#                 config_path = op.relpath(config, derivs_dir)
#                 jobs.append((derivs_dir, config_path))

#     if cores > 1:
#         with multiprocessing.Pool(cores) as p: 
#             p.starmap(__run_oxasl_worker, jobs)
#     else: 
#         for job in jobs: 
#             __run_oxasl_worker(*job)




