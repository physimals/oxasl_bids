#! /usr/bin/env python3

import argparse
import logging
import os
import sys

from .oxasl_bids import oxasl_config_from_bids, get_command_line

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
    parser.add_argument('--pipeline', choices=["oxasl", "oxford_asl"], default="oxford_asl")
    parser.add_argument('--run-fslanat', action='store_true', help="Include commands to run FSL_ANAT on structural images")
    parser.add_argument('--debug', action='store_true', help="Enable debug logging")

    args, remainder = parser.parse_known_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARN)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)

    fsl_anat_calls = []
    for config in oxasl_config_from_bids(args.bidsdir):
        if args.run_fslanat and "struct" in config:
            struct_data = config.pop("struct")
            config["fslanat"] = os.path.basename(struct_data[:struct_data.index(".nii")]) + ".anat"
            if struct_data not in fsl_anat_calls:
                # Only run fsl_anat once for each structural file
                print("fsl_anat -i %s\n" % struct_data)
                fsl_anat_calls.append(struct_data)
        print(get_command_line(config, prog=args.pipeline, extra_args=remainder) + "\n")

if __name__ == "__main__":
    main()
