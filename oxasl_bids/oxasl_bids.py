import os.path as op 
import json
import glob 
import re 
import functools
import os 
import argparse
import copy 
import warnings
import sys 
import copy 
import functools
import subprocess
import multiprocessing as mp
import shutil

import nibabel
import numpy as np
import bids
#import toblerone
#from scipy.ndimage import affine_transform

from . import utils
from .mappings import get_oxasl_config_from_metadata

__BIDS_ASL_VERSION__ = '0.0.1'

def sidecar_json(asl_dir, asl_file):
    stub = asl_file[:asl_file.index('.nii')]
    jsfile = op.join(asl_dir, stub + '.json')
    js_dict = json.load(open(jsfile, 'r'))
    return js_dict


def configuration_dir(asl_dir, asl_file):
    ddir = oxasl_dir(asl_dir, asl_file)
    cdir = op.join(ddir, 'config')
    os.makedirs(cdir, exist_ok=True)
    return cdir 


def oxasl_dir(asl_dir, asl_file):
    bids_root = asl_dir[:asl_dir.index("sub-")]
    bids_subdir = asl_dir[asl_dir.index("sub-"):]
    oxdir = op.join(bids_root, 'derivatives', "oxasl", bids_subdir)
    os.makedirs(oxdir, exist_ok=True)
    return oxdir


def derivatives_dir(asl_dir):
    return op.join(asl_dir.replace('sourcedata', 'derivatives'))


def get_asl_config(asl_file):
    """
    Identify the location of the ASL and calibration data
    and the ASL data order
    """
    options = {"asldata" : asl_file.path}
    options.update(get_oxasl_config_from_metadata(asl_file.get_metadata(), "asl"))
    metadata = asl_file.get_metadata()

    # Get the ASL context and interpret it
    ctx_filename = None
    for bids_file in asl_file.get_associations():
        if bids_file.suffix == "aslcontext":
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
        print(f"Extracting M0 from {asl_file.filename}")
        options["calib"] = asl_file.path
        options["calib_volumes"] = calib_frames
        options["asldata_volumes"] = asl_frames
        options.update(get_oxasl_config_from_metadata(asl_file.get_metadata(), "calib"))
    else: 
        # No sign of m0scan volumes in ASL context - check M0 type is separate
        # and look for it in associated files
        if metadata.get('M0Type', 'separate').lower() != "separate":
            raise utils.IncompatabilityError("JSON M0Type field not set to 'separate', but m0scan not in aslcontext")
        for bids_file in asl_file.get_associations():
            if bids_file.suffix == "m0scan":
                options["calib"] = bids_file.path
                options.update(get_oxasl_config_from_metadata(bids_file.get_metadata(), "calib"))
    return options

def get_struct_config(struct_file):
    return {"struct" : struct_file.path}

def get_calib_config(calib_file):
    # FIXME
    return {}

def json_path(niifile):
    dir, fname = op.split(niifile)
    stub = fname[:fname.index('.nii')]
    return op.join(dir, stub + '.json')

def get_oxasl_config(asl_file, sess_files):
    """
    Build configuration options from JSON metadata
    and ASL context file
    """
    options = get_asl_config(asl_file)
    for bids_file in asl_file.get_associations():
        if bids_file.suffix == "m0scan":
            # Process associated M0 file (calib or cblip)
            pass

    for struct_file in sess_files["T1w"]:
        options.update(get_struct_config(struct_file))

    for calib_file in sess_files["m0scan"]:
        options.update(get_calib_config(calib_file))

    #options.update(map_keys(json_path(options["asldata"][0]), "asl"))
    #options.update(map_keys(json_path(options["calib"][0]), "calib"))

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

def dump_to_string(arg_dict, outputdir=""):
    txt = 'oxasl '
    for key,val in arg_dict.items(): 
        if isinstance(val, list):
            val = ",".join([str(v) for v in val])

        if isinstance(val, tuple):
            if val[1] is None:
                val = op.abspath(val[0])

        if isinstance(val, bool):
            if val: 
                txt += f"--{key} "
        elif val is not None: 
            #if isinstance(val, str) and op.exists(val):
            #    val = op.relpath(val, outputdir)
            txt += f"--{key}={val} "
    return txt 

def guarded_merge(common, specific):
    merged = copy.deepcopy(specific)
    for aidx, arg_string in enumerate(common):
        arg_string = arg_string.strip('-')
        if arg_string.count('='):
            key = arg_string[:arg_string.index('=')]
            val = arg_string[arg_string.index('=')+1:]

            if key in specific:
                print(f"Warning: overwriting --{key} with value set in --common_args " + 
                f"({val} replaces {specific[key]}).")
            merged[key] = val 

        else: 
            merged[arg_string] = True

    return merged


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


def wipe_dir(path):
    shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)

def get_asl_files(bids_root):
    dataset = bids.BIDSLayout(bids_root)
    data_files = {}
    for subj in dataset.get_subjects():
        data_files[subj] = {}
        sessions = dataset.get_sessions(subject=subj)
        if not sessions:
            sessions = [None]
        for sess in sessions:
            print("Getting ASL data for subject %s in session %s" % (subj, sess))
            data_files[subj][sess] = {}
            for suffix in ("asl", "m0scan", "T1w"):
                data_files[subj][sess][suffix] = []
                for fobj in dataset.get(subject=subj, session=sess, suffix=suffix):
                    if isinstance(fobj, bids.layout.models.BIDSImageFile):
                        print("Found %s image: %s %s" % (suffix.upper(), fobj.filename, fobj.suffix))
                        data_files[subj][sess][suffix].append(fobj)

    return data_files
       
    sys.exit(1)
    asl_files = []
    for asl_dir in walk_modality_dirs(bids_root, 'perf'):
        asl_niftis = glob.glob(op.join(asl_dir, '*_asl.nii.gz'))
        asl_niftis = [ op.split(p)[1] for p in asl_niftis ]

        for asl_file in asl_niftis:
            asl_files.append((asl_dir, asl_file))
    return asl_files

def get_oxasl_configs(bids_root):
    configs = []
    asl_files = get_asl_files(bids_root)
    for subj, sessions in asl_files.items():
        for sess, sess_files in sessions.items():
            for asl_file in sess_files["asl"]:
                oxasl_args = get_oxasl_config(asl_file, sess_files)
                #oxasl_args = guarded_merge(common_args, oxasl_args)
                print("OXASL config for file %s" % asl_file.filename)
                print(dump_to_string(oxasl_args))
                configs.append(oxasl_args)
    return configs

def prepare_config_files(argv): 
    parser = argparse.ArgumentParser()
    parser.add_argument('--bidsdir', required=True, type=str)
    parser.add_argument('--common_args', required=False, 
        type=str, nargs=argparse.REMAINDER)
    parser.add_argument('--align', type=str)
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--fsl_anat', action='store_true')

    args = dict(vars(parser.parse_args(argv)))
    bids_root = args.pop('bidsdir')
    align_spc = args.pop('align')
    fsl_anat = args.pop('fsl_anat')
    overwrite = args.pop('overwrite')
    common_args = args['common_args']

    if align_spc and (align_spc != 'anat'):
        raise RuntimeError("Not implemented yet")
    if align_spc and (not fsl_anat):
        raise RuntimeError("--align must be used with --fsl_anat") 

            # oxdir = oxasl_dir(asl_dir, asl_file)
            # if overwrite: 
            #     wipe_dir(oxdir)
            # elif os.listdir(oxdir):
            #     raise RuntimeError(f"Oxasl output directory {oxdir} is not empty. Use --overwrite option.")
            # config_dir = configuration_dir(asl_dir, asl_file)
            # outpath =  op.join(config_dir, 'oxasl_config.txt')
            # rel_path_root = op.abspath(op.join(oxdir, '..'))
            # outstring = dump_to_string(oxasl_args, rel_path_root)
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


def walk_modality_dirs(bids_root, modality, datatype_dir=None):
    if datatype_dir:
        bids_root = op.join(bids_root, datatype_dir)
    all_dirs = glob.glob(op.join(bids_root, '*'))
    fltr = re.compile('sub-\d*')
    subdirs = [ p for p in all_dirs if fltr.match(op.split(p)[1]) ]
    if not subdirs: 
        raise RuntimeError(f"Did not find any subjects in {bids_root}")
    for sub in subdirs:
        for asl_dir, subdirs, files in os.walk(op.join(sub, modality)):
            dname = op.split(asl_dir)[1]
            # Exclude oxasl configuration directories 
            fltr = re.compile('sub-\d*.*_oxasl')
            if (files or subdirs) and not fltr.match(dname): 
                yield asl_dir


def execute_job(exec_dir, config_file):
    os.chdir(exec_dir)
    cmd = open(config_file, 'r').read()
    subprocess.run(cmd, shell=True)


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

def __run_oxasl_worker(at_dir, cmd_path):
    os.chdir(at_dir)
    cmd = open(cmd_path, 'r').read()
    return subprocess.run(cmd, shell=True)


def run_config_files(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('--bidsdir', required=True, type=str)
    parser.add_argument('--cores', default=1, type=int)
    args = dict(vars(parser.parse_args(argv)))
    bids_root = args['bidsdir']
    cores = args['cores']

    jobs = []  
    for derivs_dir in walk_modality_dirs(bids_root, 'perf', datatype_dir="derivatives"):
        # Look for the oxasl_directory
        all_dirs = glob.glob(op.join(derivs_dir, 'sub-*_oxasl'))
        fltr = re.compile('sub-\d*.*_oxasl')
        oxasl_dirs = [ p for p in all_dirs if fltr.match(op.split(p)[1]) ]
        for oxdir in oxasl_dirs:
            config = op.join(oxdir, 'config', 'oxasl_config.txt')
            if not op.exists(config):
                print(f"Warning: expected to find a configuration file ({config}) in {oxdir}.")
                continue
            else: 
                config_path = op.relpath(config, derivs_dir)
                jobs.append((derivs_dir, config_path))

    if cores > 1:
        with multiprocessing.Pool(cores) as p: 
            p.starmap(__run_oxasl_worker, jobs)
    else: 
        for job in jobs: 
            __run_oxasl_worker(*job)




