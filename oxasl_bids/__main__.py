#! /usr/bin/env python3

import sys 

from . import prepare_config_files, run_config_files, run_fsl_anat

help_string = ("""oxasl_bids: BIDS interface for the oxford_asl pipeline

Usage:
oxasl_bids prepare --bidsdir <bids_directory> [--common_args]
oxasl_bids run --bidsdir <bids_directory> 

"Prepare" creates an oxford_asl configuration file for each acquisition within the BIDS directory. 
Common arguments will be added in to each configuration, overriding BIDS-derived parameters.
After having used "prepare", the "run" command will execute the oxford_asl commands 
Run either command without arguments for more information.""")

def main():
    argv = sys.argv[1:]
    cmd_dict = {'prepare': prepare_config_files, 
                'run': run_config_files, 
                'fsl_anat': run_fsl_anat }
    if (len(argv)) and (argv[0] in cmd_dict): 
        cmd_dict[argv[0]](argv[1:])
    else: 
        print(help_string)

if __name__ == "__main__":
    main()
