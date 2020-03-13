import os.path as op 
import json 
import glob

def interpret_labeling_type(js_dict):
    if (js_dict['LabelingType'] != 'PASL'): 
        return True 
    else: 
        return None 


def interpret_bolus(js_dict):
    if 'PASL' in js_dict:
        return js_dict['BolusCutOffTimingSequence']
    elif ('PCASL' in js_dict) or ('CASL' in js_dict):
        return js_dict['LabelingDuration']


def calc_slicedt(js_dict):
    times = js_dict['SliceTiming']
    dts = [ pair[0]-pair[1] for pair in zip(times[1:], times[:-1]) ] 
    return sum(dts) / len(dts)


def calc_effective_echo_time(js_dict):
    # Only run this if effective echo not given 
    pass

def interpret_pedir(js_dict):
    pedir = js_dict['PhaseEncodingDirection']
    if pedir.count('-'):
        pedir = '-' + pedir.strip('-')
    return pedir

def interpret_pedir(js_dict):
    pedir = js_dict['PhaseEncodingDirection']
    if pedir.count('-'):
        pedir = '-' + pedir.strip('-')
    return pedir

in_keys = [

    # ASL options
    ('order', None),
    ('tis', 'InitialPostLabelDelay'),
    ('tes', 'EchoTime'),
    ('plds', 'InitialPostLabelDelay'),
    ('ntis', None),
    ('nplds', None),
    ('rpts', None),
    ('nphases', None),
    ('nenc', None),
    ('casl', interpret_labeling_type),
    ('bolus', interpret_bolus),
    ('slicedt', calc_slicedt),
    ('sliceband', None),
    ('artsupp', 'VascularCrushing'),

    # Structural options
    ('struc', None), 
    ('fsl_anat', None), 

    # Calibration options 
    ('calib-alpha', None),
    ('tr', None),
    ('te', 'EchoTime'),

    # Distortion correction options
    ('echospacing', 'EffectiveEchoSpacing'),
    ('pedir', interpret_pedir)

]

def map_keys(js_path):
    asldir = op.split(js_path)[0]
    js_dict = json.load(open(js_path, 'r'))
    oxasl_dict = {} 

    for (key, mapper) in in_keys:
        val = None 
        if type(mapper) is str:
            val = js_dict.get(mapper)
        elif callable(mapper) and (key in js_dict):
            val = mapper(js_dict)

        oxasl_dict[key] = val

    return oxasl_dict

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