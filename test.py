from oxasl_bids.__main__ import main
import sys 

if __name__ == "__main__":
    
    # BIDS_DIR = '/mnt/hgfs/Data/maastricht_bids'
    BIDS_DIR = 'BIDS_dataset'

# --align anat
    # common_arg_string = '--mc --bat=2 --wp --debug '
    # cmd = f'prepare --bidsdir {BIDS_DIR} --common_args ' + common_arg_string

    cmd = f'run --bidsdir {BIDS_DIR}'

    main(cmd.split())