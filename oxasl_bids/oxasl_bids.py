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
#import toblerone
#from scipy.ndimage import affine_transform

from . import utils
from .mappings import map_keys

__BIDS_ASL_VERSION__ = '0.0.1'

def sidecar_json(asl_dir, asl_file):
    stub = asl_file[:asl_file.index('.nii')]
    jsfile = op.join(asl_dir, stub + '.json')
    js_dict = json.load(open(jsfile, 'r'))
    return js_dict


def configuration_dir(asl_dir, asl_file):
    ddir = oxasl_dir(asl_dir, asl_file)
    cdir = op.join(ddir, 'bids_init')
    os.makedirs(cdir, exist_ok=True)
    return cdir 


def oxasl_dir(asl_dir, asl_file):
    stub = asl_file[:asl_file.index('asl.nii')]
    bids_root = asl_dir[:asl_dir.index("sub-")]
    bids_subdir = asl_dir[asl_dir.index("sub-"):]
    oxdir = op.join(bids_root, 'derivatives', bids_subdir, stub + 'oxasl')
    os.makedirs(oxdir, exist_ok=True)
    return oxdir


def derivatives_dir(asl_dir):
    return op.join(asl_dir.replace('sourcedata', 'derivatives'))


def extract_asl_and_calib(asl_dir, asl_file):
    stub = asl_file[:asl_file.index('asl.nii')]
    ctx_file = op.join(asl_dir, stub + 'aslcontext.tsv')
    with open(ctx_file) as f:
        ctx = [l.lower().strip() for l in f.readlines()]
    if ctx[0] != 'volume_type':
        raise utils.IncompatabilityError("First line in ASL context is not volume_type")
    ctx = ctx[1:]
    if 'cbf' in ctx:
        raise utils.IncompatabilityError("ASL context contains CBF images")
    if ('deltam' in ctx) and (('label' in ctx) or ('control' in ctx)):
        raise utils.IncompatabilityError("ASL sequence is mixed deltam and control/label")

    keys = {'control' : 0,
            'control(perev)': 0,
            'label': 1,
            'deltam': 2,
            'm0scan': 3}

    frame_codes = np.array([ keys[frame] for frame in ctx ])
    asl_frames = np.array([ i for i,code in enumerate(frame_codes) if code < 3 ])
    calib_frames = np.array([ i for i,code in enumerate(frame_codes) if code == 3 ])
    deriv_dir = derivatives_dir(asl_dir)
    js_dict = sidecar_json(asl_dir, asl_file)

    if frame_codes[asl_frames[0]] == 2:
        asl_format = 'diff'
    elif frame_codes[asl_frames[0]] > frame_codes[asl_frames[1]]:
        asl_format = 'tc'
    elif frame_codes[asl_frames[0]] < frame_codes[asl_frames[1]]:
        asl_format = 'ct'
    else: 
        raise RuntimeError("Could not interpret label-control order")

    if calib_frames.size:
        # Extract and save the ASL and calib frames separately 
        if js_dict.get('M0Type', '').lower() != "included":
            raise utils.IncompatabilityError("JSON M0Type field not set to 'included', but m0scan in aslcontext")
        asl_path = op.join(op.join(asl_dir, asl_file))
        print(f"Extracting M0 from {asl_path}")
        config_dir = configuration_dir(asl_dir, asl_file)
        img = nibabel.load(asl_path)
        series = img.get_fdata()
        asl_path = op.join(config_dir, stub + 'asl_no_M0.nii.gz')
        calib_path = op.join(config_dir, stub + 'M0.nii.gz')
        nii_class = type(img)
        asl_img = nii_class(series[...,asl_frames], img.affine,
            header=img.header)
        calib_img = nii_class(series[...,calib_frames], img.affine,
            header=img.header)
        nibabel.save(asl_img, asl_path)
        nibabel.save(calib_img, calib_path)

    else: 
        js_dict = sidecar_json(asl_dir, asl_file)
        if js_dict.get('M0Type', '').lower() != "separate":
            raise utils.IncompatabilityError("JSON M0Type field not set to 'separate', but m0scan not in aslcontext")
        calib_path = op.join(asl_dir, asl_file.replace("asl.nii", "m0scan.nii"))
        if not op.exists(calib_path):
            raise FileNotFoundError(f"Could not find M0: {calib_path}")            
        asl_path = op.join(asl_dir, asl_file)

    assert op.exists(asl_path)
    assert op.exists(calib_path)
    asl_path = op.relpath(asl_path, deriv_dir)
    calib_path = op.relpath(calib_path, deriv_dir)
    return asl_path, asl_format, calib_path


def json_path(asl_dir, asl_file):
    stub = op.split(asl_file)[1]
    stub = stub[:stub.index('.nii')]
    return op.join(asl_dir, stub + '.json')

def make_oxasl_config(asl_dir, asl_file):
    # Look for sidecar JSON
    asl_path, asl_format, calib_path = extract_asl_and_calib(asl_dir, asl_file)
    jsfile = json_path(asl_dir, asl_file)
    oxasl_args = map_keys(jsfile)
    oxasl_args['asldata'] = asl_path
    oxasl_args['calib'] = calib_path
    oxasl_args['iaf'] = asl_format

    # Post-processing (remove mutually exclusive keys)
    if oxasl_args.get('casl'):
        oxasl_args.pop('tis')
    else: 
        oxasl_args.pop('plds')

    return oxasl_args

def dump_to_string(arg_dict, outputdir):
    txt = 'oxasl '
    for key,val in arg_dict.items(): 
        if val is True: 
            txt += f"--{key} "
        elif val is not None: 
            if isinstance(val, str) and op.exists(val):
                val = op.relpath(val, outputdir)
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



def prepare_config_files(argv): 
    parser = argparse.ArgumentParser()
    parser.add_argument('--bidsdir', required=True, type=str)
    parser.add_argument('--common_args', required=True, 
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

    file_count = 0 
    for asl_dir in walk_modality_dirs(bids_root, 'perf'):
        asl_niftis = glob.glob(op.join(asl_dir, '*_asl.nii.gz'))
        asl_niftis = [ op.split(p)[1] for p in asl_niftis ]

        for asl in asl_niftis: 
            # Write back into corresponding position in derivs directory
            try:                 
                oxasl_args = make_oxasl_config(asl_dir, asl)
                oxasl_args = guarded_merge(common_args, oxasl_args)
                oxdir = oxasl_dir(asl_dir, asl)
                if overwrite: 
                    wipe_dir(oxdir)
                elif os.listdir(oxdir):
                    raise RuntimeError(f"Oxasl output directory {oxdir} is not empty. Use --overwrite option.")
                config_dir = configuration_dir(asl_dir, asl)
                outpath =  op.join(config_dir, 'oxasl_config.txt')
                rel_path_root = op.abspath(op.join(oxdir, '..'))
                outstring = dump_to_string(oxasl_args, rel_path_root)
                outstring += "--output={} ".format(oxdir)

                if fsl_anat: 
                    anatdir = asl_dir.replace('sourcedata', 'derivatives')
                    anatdir = anatdir.replace('asl', 'anat')
                    anatdirs = glob.glob(op.join(anatdir, 'sub-*T1w.anat'))
                    if len(anatdirs) > 1:
                        raise RuntimeError(f"Found multiple fsl_anat dirs in {anatdir}")
                    elif not len(anatdirs):
                        raise RuntimeError(f"Did not find fsl_anat dir in {anatdir}")
                    fslanat = op.relpath(anatdirs[0], rel_path_root)
                    outstring += "--fslanat={} ".format(fslanat)

                if align_spc and fsl_anat: 
                    t1 = op.join(anatdirs[0], 'T1_biascorr_brain.nii.gz')
                    align_img = prepare_alignment_target(t1, op.join(asl_dir, asl))
                    align_path = op.join(config_dir, 'oxasl_bids_common_space.nii.gz')
                    nibabel.save(align_img, align_path)
                    align_mat = utils._world_to_FLIRT(t1, align_path, np.eye(4))
                    align_mat_path = op.join(config_dir, 'oxasl_bids_common_space_flirt.txt')
                    rel_mat_path, rel_algn_path = [ op.relpath(p, rel_path_root) 
                        for p in [ align_mat_path, align_path ] ] 
                    np.savetxt(align_mat_path, align_mat)
                    outstring += f' --output-custom={rel_algn_path} --output-custom-mat={rel_mat_path}'

                file_count += 1     
                outstring += ' --overwrite'
                with open(outpath, 'w') as f: 
                    f.write(outstring)

            except utils.IncompatabilityError as e:
                print(f"Warning: incompatible parameters found for {asl} - skipping")
                print(e)
                continue
            except FileNotFoundError as e: 
                print(f"Warning: missing file for {asl} - skipping")
                print(e)
                continue
            except Exception as e: 
                raise e 

    print(f"Wrote {file_count} configuration files in {op.join(bids_root, 'derivatives')}")


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
            config = op.join(oxdir, 'bids_init', 'oxasl_config.txt')
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




