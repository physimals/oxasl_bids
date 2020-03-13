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
import subprocess

import nibabel 
import bids
import numpy as np 

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
    oxdir = op.join(asl_dir.replace('sourcedata', 'derivatives'),
        stub + 'oxasl')
    os.makedirs(oxdir, exist_ok=True)
    return oxdir


def derivatives_dir(asl_dir):
    return op.join(asl_dir.replace('sourcedata', 'derivatives'))


def extract_asl_and_calib(asl_dir, asl_file):
    stub = asl_file[:asl_file.index('asl.nii')]
    ctx_file = op.join(asl_dir, stub + 'ASLContext.tsv')
    ctx = open(ctx_file).read().split()
    if 'CBF' in ctx:
        raise utils.IncompatabilityError("ASLContext contains CBF images")
    if ('DeltaM' in ctx) and (('Label' in ctx) or ('Control' in ctx)):
        raise utils.IncompatabilityError("ASL sequence is mixed DeltaM and Control/Label")

    keys = {'Control' : 0,
            'Control(PErev)': 0,
            'Label': 1, 
            'DeltaM': 2,
            'M0Scan': 3}

    indices = np.array([ keys[frame] for frame in ctx ])
    asl_frames = np.array([ i for i,idx in enumerate(indices) if idx < 3 ])
    calib_frames = np.array([ i for i,idx in enumerate(indices) if idx == 3 ])
    deriv_dir = derivatives_dir(asl_dir)
    js_dict = sidecar_json(asl_dir, asl_file)

    if indices[asl_frames[0]] == 3:
        asl_format = 'diff'
    elif indices[asl_frames[0]] > indices[asl_frames[1]]:
        asl_format = 'tc'
    elif indices[asl_frames[0]] < indices[asl_frames[1]]:
        asl_format = 'ct'
    else: 
        raise RuntimeError("Could not interpret label-control order")

    if calib_frames: 
        # Extract and save the ASL and calib frames separately 
        assert (js_dict['M0'] is True), 'JSON M0 field does not match ASLContext' 
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
        if not isinstance(js_dict['M0'], str): 
            raise utils.IncompatabilityError("M0 field in JSON must point to separate acquisition")
        calib_path = op.join(asl_dir, js_dict['M0'])
        if not op.exists(calib_path):
            raise FileNotFoundError(f"Could not find M0: {calib_path}")            
        asl_path = asl_file

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
    oxasl_args['iaf']

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
    
    vol = np.zeros(dims, dtype=np.int8)
    align_img = nibabel.Nifti2Image(vol, affine=new_affine)

    return align_img



def prepare_config_files(argv): 
    parser = argparse.ArgumentParser()
    parser.add_argument('--bidsdir', required=True, type=str)
    parser.add_argument('--common_args', required=True, 
        type=str, nargs=argparse.REMAINDER)
    parser.add_argument('--align', type=str)

    args = dict(vars(parser.parse_args(argv)))
    bids_root = args.pop('bidsdir')
    align_spc = args.pop('align')
    common_args = args['common_args']

    file_count = 0 
    for asl_dir in walk_asl_dirs(bids_root, 'sourcedata'):
        asl_niftis = glob.glob(op.join(asl_dir, '*_asl.nii.gz'))
        asl_niftis = [ op.split(p)[1] for p in asl_niftis ]

        # Grab the structural image if needed for alignment
        if align_spc: 
            align_target = get_alignment_path(align_spc, asl_dir)

        for asl in asl_niftis: 
            # Write back into corresponding position in derivs directory
            try:                 
                oxasl_args = make_oxasl_config(asl_dir, asl)
                oxasl_args = guarded_merge(common_args, oxasl_args)
                oxdir = oxasl_dir(asl_dir, asl)
                config_dir = configuration_dir(asl_dir, asl)
                outpath =  op.join(config_dir, 'oxasl_config.txt')
                rel_path_root = op.abspath(op.join(oxdir, '..'))
                outstring = dump_to_string(oxasl_args, rel_path_root)

                if align_spc: 
                    align_img = prepare_alignment_target(align_target, asl)
                    align_path = op.join(oxdir, 'oxasl_bids_common_space.nii.gz')
                    nibabel.save(align_img, align_path)
                    align_mat = utils._world_to_FLIRT(align_target, align_path, np.eye(4))
                    align_mat_path = op.join(oxdir, 'oxasl_bids_common_space_flirt.txt')
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
            except FileNotFoundError: 
                print(f"Warning: could not find JSON for {asl} - skipping")
                continue
            except Exception as e: 
                raise e 

    print(f"Wrote {file_count} configuration files in {op.join(bids_root, 'derivatives')}")


def walk_asl_dirs(bids_root, datatype):
    all_dirs = glob.glob(op.join(bids_root, datatype, '*'))
    subdirs = [ p for p in all_dirs if re.match('sub-\d*', op.split(p)[1]) ]
    for sub in subdirs:
        for asl_dir, subdirs, files in os.walk(op.join(sub, 'asl')):
            dname = op.split(asl_dir)[1]
            # Exclude oxasl configuration directories 
            filter = 'sub-\d*.*_oxasl'
            if (files or subdirs) and not re.match(filter, dname): 
                yield asl_dir

def execute_job(exec_dir, config_file):
    os.chdir(exec_dir)
    cmd = open(config_file, 'r').read()
    subprocess.run(cmd, shell=True)



def run_config_files(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('--bidsdir', required=True, type=str)
    args = dict(vars(parser.parse_args(argv)))

    bids_root = args['bidsdir']

    jobs = []  
    for asl_dir in walk_asl_dirs(bids_root, 'derivatives'): 
        # Look for the oxasl_directory
        all_dirs = glob.glob(op.join(asl_dir, 'sub-*_oxasl'))
        filter = 'sub-\d*.*_oxasl'
        oxasl_dirs = [ p for p in all_dirs if re.match(filter, op.split(p)[1]) ]
        for oxdir in oxasl_dirs:
            config = op.join(oxdir, 'bids_init', 'oxasl_config.txt')
            if not op.exists(config):
                print(f"Warning: expected to find a configuration file ({config}) in {oxdir}.")
                continue
            else: 
                config = op.relpath(config, asl_dir)
                jobs.append((asl_dir, config))

    for job in jobs: 
        execute_job(*job)








    # x = asl 
    # key_vals = []
    # while len(x):  
    #     ks = 0  
    #     ke = x.index('-') 
    #     key = x[ks:ke] 
    #     vs = ke + 1 
    #     ve = x.index('_') 
    #     val = x[vs:ve] 
    #     key_vals.append((key,val))
    #     if not x[ve+1:].count('-'): 
    #         break 
    #     else: 
    #         x = x[ve+1:]