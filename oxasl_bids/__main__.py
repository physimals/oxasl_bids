"""
oxasl_bids: BIDS interface for the oxford_asl pipeline

Usage:
oxasl_bids script --bidsdir <bids_directory> [... additional oxford_asl arguments]
oxasl_bids args --img <img file> --img-type <type>
oxasl_bids bidsout --oxasl-output <output folder>

'script' outputs commands to run oxasl or oxford_asl on a BIDS data set
'args' translates metadata from an image file with associated JSON sidecar into oxford_asl command line arguments
'bidsout' converts an oxford_asl output directory to BIDS format
"""

import argparse
import json
import logging
import os
import sys

import nibabel as nib

from .oxasl_bids import oxford_asl_to_bids, oxasl_config_from_bids, get_command_line, get_fslanat_command, get_output_as_bids_command
from .mappings import get_oxasl_config_from_metadata

def _parse_args(args):
    # FIXME not used at present
    ret = {}
    skip = False
    for idx, arg in enumerate(args):
        if skip:
            skip = False
            continue

        if arg.startswith("-"):
            arg = arg.strip("-")
            parts = arg.split("=", 1)
            if len(parts) == 2:
                k, v = parts
            elif idx < len(args)-1 and not args[idx+1].startswith("-"):
                k = parts[0]
                v = args[idx+1]
                skip = True
            else:
                k = parts[0]
                v = ""
            ret[k.strip()] = v.strip()
        else:
            raise RuntimeError("Invalid argument: %s (does not start with -)" % arg)
    return ret

def main():
    parser = argparse.ArgumentParser()
    group = parser.add_argument_group('General options')
    group.add_argument('command', help="Command to run", choices=["script", "args", "bidsout"])
    group.add_argument('--pipeline', help="Target pipeline", choices=["oxasl", "oxford_asl"], default="oxford_asl")
    group.add_argument('--debug', help="Enable debug logging", action='store_true')
    group = parser.add_argument_group('Script mode options')
    group.add_argument('--bidsdir', help="Path to BIDS data set")
    group.add_argument('--run-fslanat', help="Include commands to run FSL_ANAT on structural images", action='store_true')
    group = parser.add_argument_group('Args mode options')
    group.add_argument('--img', help="Path to image file with matching JSON sidecar")
    group.add_argument('--img-type', help="Image type", choices=["asl", "calib", "cblip"])
    group = parser.add_argument_group('Bids output mode options (--bidsdir required)')
    group.add_argument('--oxasl-output', help="Path to oxford ASL or oxasl output")
    group.add_argument('--bids-output', help="Directory to store BIDS output in (if not using --merge-source)")
    group.add_argument('--merge-source', help="Merge output with source BIDS dataset (specified with --bidsdir)", action='store_true')
    group.add_argument('--subject', help="Subject ID for output")
    group.add_argument('--session', help="Session ID for output")

    args, remainder = parser.parse_known_args()
    custom_args = _parse_args(remainder)

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
    elif args.command == "bidsout":
        do_bidsout(args, remainder)

def do_script(args, remainder):
    """
    Generate a script to run oxford_asl or oxasl on a BIDS data set

    The script will include an oxasl/oxford_asl call for each instance
    of ASL data found, in addition there may be fsl_anat calls for
    processing structural data
    """
    fsl_anat_calls = []
    for config in oxasl_config_from_bids(args.bidsdir):
        options = config["options"]
        if args.run_fslanat and "struct" in options:
            fslanat_cmd = get_fslanat_command(options)
            if fslanat_cmd not in fsl_anat_calls:
                print(fslanat_cmd)
                fsl_anat_calls.append(fslanat_cmd)
        print(get_command_line(options, prog=args.pipeline, extra_args=remainder) + "\n")
        if args.bids_output or args.merge_source:
            output_asl_bids_cmd = get_output_as_bids_command(args, config)
            print(output_asl_bids_cmd)

def do_args(args, remainder):
    """
    Get the oxasl/oxford_asl configuration options for a BIDS dataset
    """
    img_shape = nib.load(args.img).shape
    json_filename = args.img[:args.img.index(".nii")] + ".json"
    with open(json_filename, "r") as f:
        json_metadata = json.load(f)

    options = {args.img_type : args.img}
    options.update(get_oxasl_config_from_metadata(json_metadata, args.img_type, img_shape))
    print(get_command_line(options, prog="", extra_args=remainder))

def do_bidsout(args, remainder):
    """
    Take an existing oxford_asl/oxasl output and convert it to BIDS

    The results may either be written as an independent BIDS data set
    or merged with the source BIDS data set as a deriviative
    """
    if not args.bids_output and not args.merge_source:
        raise RuntimeError("If --bids-output is not provided --merge-source must be specified")
    if  args.bids_output and args.merge_source:
        raise RuntimeError("Can't specify --bids-output and --merge-source at the same time")
    print(f"Generating BIDS output for oxford_asl output at {args.oxasl_output}")
    oxford_asl_to_bids(args.oxasl_output, args.bidsdir, subject=args.subject, session=args.session, bids_output_dir=args.bids_output)

if __name__ == "__main__":
    main()
