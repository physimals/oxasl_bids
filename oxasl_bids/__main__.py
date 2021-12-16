#! /usr/bin/env python3

import argparse
import json
import logging
import os
import sys

import nibabel as nib

from .oxasl_bids import oxasl_config_from_bids, get_command_line
from .mappings import get_oxasl_config_from_metadata

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
    parser.add_argument('command', help="Command to run", choices=["script", "args"])
    parser.add_argument('--pipeline', help="Target pipeline", choices=["oxasl", "oxford_asl"], default="oxford_asl")
    parser.add_argument('--debug', help="Enable debug logging", action='store_true')
    group = parser.add_argument_group('Script mode options')
    group.add_argument('--bidsdir', help="Path to BIDS data set")
    group.add_argument('--run-fslanat', help="Include commands to run FSL_ANAT on structural images", action='store_true')
    group = parser.add_argument_group('Args mode options')
    group.add_argument('--img', help="Path to image file with matching JSON sidecar")
    group.add_argument('--img-type', help="Image type", choices=["asl", "calib", "cblip"])
    
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

    if args.command == "script":
        do_script(args, remainder)
    elif args.command == "args":
        do_args(args, remainder)

def do_script(args, remainder):
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

def do_args(args, remainder):
    img_shape = nib.load(args.img).shape
    json_filename = args.img[:args.img.index(".nii")] + ".json"
    with open(json_filename, "r") as f:
        json_metadata = json.load(f)

    options = get_oxasl_config_from_metadata(json_metadata, args.img_type, img_shape)
    print(get_command_line(options, prog="", extra_args=remainder))

if __name__ == "__main__":
    main()
