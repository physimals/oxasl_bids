"""
oxasl_bids: BIDS interface for the oxasl pipeline

Usage:
oxasl_bids bids --bidsdir <bids_directory> [... additional oxasl arguments]
oxasl_bids img --img <img file> --img-type <type>
oxasl_bids bidsout --oxasl-output <output folder>

'bids' outputs commands to run oxasl on a BIDS data set
'img' output a command to run oxasl on a single ASL image with JSON metadata
'bidsout' converts an oxasl output directory to BIDS format
"""

import argparse
import logging
import sys

from .oxasl_bids import oxasl_output_to_bids, oxasl_config_from_bids, get_fslanat_command, get_output_as_bids_command, get_oxasl_command_line
from .mappings import oxasl_config_from_metadata
from . import utils

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

class ArgumentParser(argparse.ArgumentParser):
    def __init__(self, **kwargs):
        argparse.ArgumentParser.__init__(self, prog="oxasl_bids", add_help=True, **kwargs)
        group = self.add_argument_group('General options')
        group.add_argument('mode', help="Run mode", choices=["bids", "img", "bidsout"])
        group.add_argument('--debug', help="Enable debug logging", action='store_true')
        group = self.add_argument_group('Script mode options')
        group.add_argument('--bidsdir', help="Path to BIDS data set")
        group.add_argument('--run-fslanat', help="Include commands to run FSL_ANAT on structural images", action='store_true')
        group = self.add_argument_group('Args mode options')
        group.add_argument('--img', help="Path to image file with matching JSON sidecar")
        group.add_argument('--img-type', help="Image type", choices=["asl", "calib", "cblip", "struc"])
        group = self.add_argument_group('Bids output mode options (--bidsdir required)')
        group.add_argument('--oxasl-output', help="Path to oxasl output")
        group.add_argument('--bids-output', help="Directory to store BIDS output in (if not using --merge-source)")
        group.add_argument('--merge-source', help="Merge output with source BIDS dataset (specified with --bidsdir)", action='store_true')
        group.add_argument('--subject', help="Subject ID for output")
        group.add_argument('--session', help="Session ID for output")

def main():
    parser = ArgumentParser()
    args, remainder = parser.parse_known_args()
    custom_args = _parse_args(remainder)

    setup_logging(args)
    
    if args.mode == "bids":
        do_bids(args, remainder)
    elif args.mode == "img":
        do_img(args, remainder)
    elif args.mode == "bidsout":
        do_bidsout(args, remainder)

def setup_logging(args):
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.WARN)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)

def do_bids(args, remainder):
    """
    Generate a script to run oxasl on all sessions in a BIDS data set

    The script will include an oxasl call for each instance
    of ASL data found, in addition there may be fsl_anat calls for
    processing structural data
    """
    fsl_anat_calls = []
    for config in oxasl_config_from_bids(args.bidsdir):
        options = config["options"]
        if args.run_fslanat and "struc" in options:
            fslanat_cmd = get_fslanat_command(options)
            if fslanat_cmd not in fsl_anat_calls:
                print(fslanat_cmd)
                fsl_anat_calls.append(fslanat_cmd)
        print(get_oxasl_command_line(options, extra_args=remainder) + "\n")
        if args.bids_output or args.merge_source:
            output_asl_bids_cmd = get_output_as_bids_command(args, config)
            print(output_asl_bids_cmd)

def do_img(args, remainder):
    """
    Get the oxasl command for an ASL image file with JSON metadata

    The command will not contain any options for structural or calibration data
    """
    nii, metadata = utils.load_img(args.img)
    options = {args.img_type : args.img}
    options.update(oxasl_config_from_metadata(metadata, args.img_type))
    print(get_oxasl_command_line(options, extra_args=remainder, prog=""))

def do_bidsout(args, remainder):
    """
    Take an existing oxasl output and convert it to BIDS

    The results may either be written as an independent BIDS data set
    or merged with the source BIDS data set as a deriviative
    """
    if not args.bids_output and not args.merge_source:
        raise RuntimeError("If --bids-output is not provided --merge-source must be specified")
    if  args.bids_output and args.merge_source:
        raise RuntimeError("Can't specify --bids-output and --merge-source at the same time")
    print(f"Generating BIDS output for oxasl output at {args.oxasl_output}")
    oxasl_output_to_bids(args.oxasl_output, args.bidsdir, subject=args.subject, session=args.session, bids_output_dir=args.bids_output)

if __name__ == "__main__":
    main()
