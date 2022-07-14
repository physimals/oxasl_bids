"""
Runs the converter on the ASL examples in the BIDS examples
Github module
"""
import os

from oxasl_bids.oxasl_bids import oxasl_config_from_bids

EXAMPLES_DIR = "/home/martin/code/thirdparty/bids-examples"

EXAMPLES = {
    "asl001" : {
        # pre-subtracted PCASL single PLD with embedded M0 in first volume
        "asldata" : ("perf/sub-Sub103_asl.nii.gz", 1),
        "iaf" : "diff",
        "tes" : 0.010528,
        "casl" : True,
        "bolus" : 1.450,
        "plds" : 2.025,
        "calib" : ("perf/sub-Sub103_asl.nii.gz", 0),
        "tr" : 4.886,
        "te" : 0.010528,
        "struct" : ("anat/sub-Sub103_T1w.nii.gz", None),
    },
    "asl002" : {
        # PCASL single PLD CT pairs 2D acquisition with separate M0
        "iaf" : "ct",
        "tes" : 0.015,
        "casl" : True,
        "bolus" : 1.8,
        "plds" : 2.0,
        "slicedt" : 0.0385,
        "asldata" : ("perf/sub-Sub103_asl.nii.gz", None),
        "calib" : ("perf/sub-Sub103_m0scan.nii.gz", None),
        "struct" : ("anat/sub-Sub103_T1w.nii.gz", None),
    },
    "asl003" : {
        # Multi-PLD PASL TC pairs with separate M0
        "iaf" : "tc",
        "tes" : 0.01192,
        "tis" : [0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3],
        "asldata" : ("perf/sub-Sub1_asl.nii.gz", None),
        "calib" : ("perf/sub-Sub1_m0scan.nii.gz", None),
        "struct" : ("anat/sub-Sub1_T1w.nii.gz", None),
    },
    "asl004" : {
        # 48 TC pairs 2D acquisition multi PLD with 8 repeats per PLD with separate M0
        # Note we do not handle repeats at present so you end up with a long list of PLDs
        "iaf" : "tc",
        "slicedt" : 0.0452,
        "bolus" : 1.4,
        "casl" : True,
        "plds" : [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
                  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                  0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75,
                  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                  1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25,
                  1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
        "asldata" : ("perf/sub-Sub1_asl.nii.gz", None),
        "calib" : ("perf/sub-Sub1_m0scan.nii.gz", None),
        "cblip" : ("fmap/sub-Sub1_dir-pa_m0scan.nii.gz", None),
        "echospacing" : 0.00095,
        "pedir" : "y",
        "struct" : ("anat/sub-Sub1_T1w.nii.gz", None),
    },
    "asl005" : {
        # 8 CT pairs single PLD PCASL separate M0
        "iaf" : "ct",
        "bolus" : 1.8,
        "plds" : 2.0,
        "asldata" : ("perf/sub-Sub103_asl.nii.gz", None),
        "calib" : ("perf/sub-Sub103_m0scan.nii.gz", None),
        "struct" : ("anat/sub-Sub103_T1w.nii.gz", None),
    },
}

def check_matches(config, expected_config):
    for key, value in expected_config.items():
        actual_value = config[key]
        if isinstance(actual_value, list) and len(actual_value) == 1:
            actual_value = actual_value[0]

        if isinstance(value, float):
            assert abs(value - actual_value) < 0.00001
        elif isinstance(value, tuple):
            assert(actual_value.endswith(value[0]))
            if value[1] is not None:
                assert(config[key + "_volumes"] == value[1])
        elif isinstance(value, list):
            assert(len(value) == len(actual_value))
            for v, av in zip(value, actual_value):
                assert abs(v - av) < 0.00001 
        else:
            assert(value == actual_value)

def check_bids_example(subdir):
    expected_config = EXAMPLES[subdir]
    bids_root = os.path.join(EXAMPLES_DIR, subdir)
    configs = oxasl_config_from_bids(bids_root)
    assert(len(configs) == 1)
    check_matches(configs[0], expected_config)

def test_asl001():
    check_bids_example("asl001")

def test_asl002():
    check_bids_example("asl002")

def test_asl003():
    check_bids_example("asl003")

def test_asl004():
    check_bids_example("asl004")

def test_asl005():
    check_bids_example("asl005")


