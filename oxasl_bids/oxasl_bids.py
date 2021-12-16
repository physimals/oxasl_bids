"""
OXASL_BIDS: Maps BIDS data sets on to oxasl/oxford_asl options
"""
import os.path as op
import os
import copy
import logging

import nibabel
import numpy as np
import bids

from . import utils
from .mappings import get_oxasl_config_from_metadata

LOG = logging.getLogger(__name__)

def get_command_line(options, prog="oxasl", extra_args=[]):
    """
    :param options: Dictionary of options derived from BIDS

    :return: String suitable to pass to oxasl or oxford_asl
    """
    txt = prog + ' '
    for key,val in options.items():
        if key == "asldata":
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

        if isinstance(val, tuple):
            if val[1] is None:
                val = op.abspath(val[0])

        if isinstance(val, bool):
            if val: 
                txt += f"{key} "
        elif val is not None: 
            #if isinstance(val, str) and op.exists(val):
            #    val = op.relpath(val, outputdir)
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
    options = {"asldata" : op.abspath(asl_file.path)}
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
            raise utils.IncompatabilityError("JSON M0Type field not set to 'included', but m0scan in aslcontext")
        LOG.debug(f"Extracting M0 from {asl_file.filename}")
        options["calib"] = op.abspath(asl_file.path)
        options["calib_volumes"] = calib_frames
        options["asldata_volumes"] = asl_frames
        options.update(get_oxasl_config_from_metadata(asl_file.get_metadata(), "calib"))
    else: 
        # No sign of m0scan volumes in ASL context - check M0 type is separate
        # and look for it in associated files
        if metadata.get('M0Type', 'separate').lower() != "separate":
            raise utils.IncompatabilityError("JSON M0Type field not set to 'separate', but m0scan not in aslcontext")
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
        ret.update(get_oxasl_config_from_metadata(m0_file.get_metadata(), "cblip"))
        if "totalreadouttime" in ret:
            # We need the effective echo spacing instead
            readouttime = ret.pop("totalreadouttime")
            pedir = ret.get("pedir", None)
            if not pedir:
                LOG.warn("Found total readout time for cblip image but no PE dir - cannot calculate effective echo spacing")
            else:
                img_dims = m0_file.get_image().shape
                if "x" in pedir:
                    size = img_dims[0]
                elif "y" in pedir:
                    size = img_dims[1]
                elif "z" in pedir:
                    size = img_dims[2]
                # FIXME what if pedir invalid
                ret["echospacing"] = readouttime / (size-1)
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

    # Find out what field we're using for timings based on whether it's PCASL or PASL
    if options.get('casl', False):
        ttype, otype = 'plds', 'tis'
    else:
        ttype, otype = 'tis', 'plds'

    # Remove unused timings field and make sure it is a list
    options.pop(otype, None)
    if isinstance(options[ttype], (int, float)):
        options[ttype] = [options[ttype]]
    if options.get('iaf', 'diff') != 'diff':
        # With TC pairs the timings will be repeated. We don't care about the order since
        # no valid ASL sequence will have different timings for tag and control.
        options[ttype] = [options[ttype][idx] for idx in range(0, len(options[ttype]), 2)]
    
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




