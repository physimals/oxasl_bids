#! /usr/bin/env python3

import argparse
import sys 

from .oxasl_bids import get_oxasl_configs

help_string = ("""oxasl_bids: BIDS interface for the oxford_asl pipeline

Usage:
oxasl_bids prepare --bidsdir <bids_directory> [--common_args]
oxasl_bids run --bidsdir <bids_directory> 

"Prepare" creates an oxford_asl configuration file for each acquisition within the BIDS directory. 
Common arguments will be added in to each configuration, overriding BIDS-derived parameters.
After having used "prepare", the "run" command will execute the oxford_asl commands 
Run either command without arguments for more information.""")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bidsdir', required=True, type=str)
    #parser.add_argument('--common_args', required=False, 
    #    type=str, nargs=argparse.REMAINDER)
    #parser.add_argument('--align', type=str)
    #parser.add_argument('--overwrite', action='store_true')
    #parser.add_argument('--fsl_anat', action='store_true')

    args = dict(vars(parser.parse_args()))
    bids_root = args.pop('bidsdir')
    get_oxasl_configs(bids_root)
    
    #align_spc = args.pop('align')
    #fsl_anat = args.pop('fsl_anat')
    #overwrite = args.pop('overwrite')
    #common_args = args['common_args']

    #if align_spc and (align_spc != 'anat'):
    #    raise RuntimeError("Not implemented yet")
    #if align_spc and (not fsl_anat):
    #    raise RuntimeError("--align must be used with --fsl_anat") 

    
if __name__ == "__main__":
    main()
