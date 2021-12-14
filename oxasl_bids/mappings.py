import os.path as op 
import json 
import glob

def is_casl(js_dict, oxasl_dict):
    label_type = js_dict.get('LabelingType', js_dict['ArterialSpinLabelingType'])
    if label_type != 'PASL': 
        return True

def calc_bolus(js_dict, oxasl_dict):
    if not oxasl_dict.get("casl", False) and "BolusCutOffTImingSequence" in js_dict:
        return js_dict['BolusCutOffTimingSequence']
    elif 'LabelingDuration' in js_dict:
        return js_dict['LabelingDuration']

def calc_slicedt(js_dict, oxasl_dict):
    if "SliceTiming" in js_dict:
        times = js_dict['SliceTiming']
        dts = [ pair[0]-pair[1] for pair in zip(times[1:], times[:-1]) ] 
        return sum(dts) / len(dts)
    else:
        return None

def calc_effective_echo_time(js_dict, oxasl_dict):
    # Only run this if effective echo not given 
    pass

def interpret_pedir(js_dict, oxasl_dict):
    pedir = js_dict['PhaseEncodingDirection']
    if pedir.count('-'):
        pedir = '-' + pedir.strip('-')
    return pedir

def interpret_pedir(js_dict, oxasl_dict):
    pedir = js_dict['PhaseEncodingDirection']
    if pedir.count('-'):
        pedir = '-' + pedir.strip('-')
    return pedir

in_keys = {

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
        ('casl', is_casl),
        ('bolus', calc_bolus),
        ('slicedt', calc_slicedt),
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
        ('echospacing', 'EffectiveEchoSpacing'),
        ('pedir', interpret_pedir),
    ],
}

def get_oxasl_config_from_metadata(json_metadata, filetype):
    """
    Get the relevant oxasl options from JSON metadata
    
    :param json_filename: Path to JSON metadata file
    :param filetype: Type of file metadata describes: asl, calib, cblip, struct
    """
    #print("mapping keys for %s" % filetype)
    oxasl_config = {} 
    #with open(json_filename, "r") as f:
    #    json_metadata = json.load(f)

    for oxasl_key, json_key in in_keys[filetype]:
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

# def interpret_m0(js_dict):
#     try: 
#         M0 = js_dict['M0']
#     except KeyError: 
#         try:
#             M0 = js_dict['MZero']
#         except KeyError: 
#             raise RuntimeError("Could not find M0 / MZero field")

#     if isinstance(M0, float) or isinstance(M0, int):
#         raise RuntimeError("Only calibration with M0 image currently supported")
#     elif isinstance(M0, str):
#         return M0
#     elif isinstance(M0, bool):
#         # FIXME: need to extract the M0 volume from within the ASL timeseries
#         raise RuntimeError("Not implemented yet")
#     else: 
#         raise RuntimeError("Could not interpret M0")