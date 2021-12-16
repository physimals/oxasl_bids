"""
OXASL_BIDS: Defines mapping of BIDS JSON metadata keys to oxasl options
"""

import logging

LOG = logging.getLogger(__name__)

def _is_casl(js_dict, oxasl_dict):
    label_type = js_dict.get('LabelingType', js_dict['ArterialSpinLabelingType'])
    if label_type != 'PASL': 
        return True

def _calc_bolus(js_dict, oxasl_dict):
    if not oxasl_dict.get("casl", False) and "BolusCutOffTImingSequence" in js_dict:
        return js_dict['BolusCutOffTimingSequence']
    elif 'LabelingDuration' in js_dict:
        return js_dict['LabelingDuration']

def _calc_slicedt(js_dict, oxasl_dict):
    if "SliceTiming" in js_dict:
        times = js_dict['SliceTiming']
        dts = [ pair[0]-pair[1] for pair in zip(times[1:], times[:-1]) ] 
        return sum(dts) / len(dts)
    else:
        return None

def _interpret_pedir(js_dict, oxasl_dict):
    dir_map = {"i" : "x", "j" : "y", "k" : "z"}
    pedir = js_dict['PhaseEncodingDirection']
    oxasl_pedir = dir_map[pedir.strip("-")]
    if pedir.count('-'):
        oxasl_pedir = '-' + oxasl_pedir
    return oxasl_pedir

OXASL_JSON_MAPPINGS = {

    "asl" : [
        # ASL options
        ('order', None),
        ('tis', 'InitialPostLabelDelay'),
        ('tis', 'PostLabelingDelay'),
        ('plds', 'InitialPostLabelDelay'),
        ('plds', 'PostLabelingDelay'),
        ('tes', 'EchoTime'),
        ('ntis', None),
        ('nplds', None),
        ('rpts', None),
        ('nphases', None),
        ('nenc', None),
        ('casl', _is_casl),
        ('bolus', _calc_bolus),
        ('slicedt', _calc_slicedt),
        ('sliceband', None),
        ('artsupp', 'VascularCrushing'),
    ],

    "struct" : [
        # Structural options
        ('struc', None), 
        ('fsl_anat', None), 
    ],

    "calib" : [
        # Calibration options 
        ('calib-alpha', None),
        ('tr', "RepetitionTimePreparation"),
        ('tr', "RepetitionTime"),
        ('te', 'EchoTime'),
    ],

    "cblip" : [
        # Distortion correction options
        ('echospacing', "EffectiveEchoSpacing"),
        ('totalreadouttime', "TotalReadoutTime"),
        ('pedir', _interpret_pedir),
    ],
}

def get_oxasl_config_from_metadata(json_metadata, filetype):
    """
    Get the relevant oxasl options from JSON metadata
    
    :param json_filename: Path to JSON metadata file
    :param filetype: Type of file metadata describes: asl, calib, cblip, struct
    """
    oxasl_config = {} 
    #with open(json_filename, "r") as f:
    #    json_metadata = json.load(f)

    for oxasl_key, json_key in OXASL_JSON_MAPPINGS[filetype]:
        #print("Looking for %s->%s" % (json_key, oxasl_key))
        val = None
        if type(json_key) is str and json_key in json_metadata:
            #print("json dict? %s" % json_metadata.get(json_key, None))
            val = json_metadata.get(json_key)
        elif callable(json_key):
            val = json_key(json_metadata, oxasl_config)

        if val is not None:
            oxasl_config[oxasl_key] = val

    return oxasl_config
