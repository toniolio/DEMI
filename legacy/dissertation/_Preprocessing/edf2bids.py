########################################
## Raw data to BIDS conversion script ##
########################################

__author__ = "Austin Hurst"

# This script is for converting the raw EEG data from the DEMI project into the
# BIDS standard, which keeps everything neatly organized and keeps all relevant
# metadata from the project in one place. It should also help with archiving
# the data for future and making other analysis scripts easier to write.



### Import required packages ###

import os
import re
import json

import mne
from mne_bids import write_raw_bids, BIDSPath
from mne_bids.utils import _write_json



### Parameters for script ###

raw_data_path = "edfs"
bids_root= "bids"
taskname = "tracelab"

raw_data_ext = ".edf"
id_num_pattern = r"demi_(\d+)[\s\.]" # ignores split files for now


### Additional conversion metadata ###

# Map of EDF annotations to MNE event codes
event_map = {
    '28': 28,
    '30': 30,
    '44': 44,
    '46': 46,
    '60': 60,
    '62': 62,
    'file start': 100
}

# Map of MNE event codes to event names
event_name_map = {
    28: "stim_on",
    30: "red_on",
    44: "trace_start",
    46: "trace_end",
    60: "accuracy_submit",
    62: "vividness_submit",
    100: "file start"
}

# Type overrides for non-EEG channels
channel_type_overrides = {
    'HEO': "eog",
    'VEO': "eog",
    'EMG-L': "emg",
    'EMG-A': "emg"
}

# Additional metadata for BIDS sidecar files
# NOTE: the 'None's here get filled in with values from the file being
# converted. They're just in the dict to specify the order of the fields.
metadata = {
    'TaskName': "tracelab",
    'TaskDescription': "[Description goes here]",
    'SamplingFrequency': None,
    'Manufacturer': "Compumedics Neuroscan",
    'ManufacturersModelName': "SynAmps 2/RT",
    'CapManufacturer': "Compumedics Neuroscan",
    'CapManufacturersModelName': "Quik-Cap 64ch [need to confirm]",
    'SoftwareVersions': "Compumedics CURRY 7",
    'EEGChannelCount': None,
    'EOGChannelCount': None,
    'ECGChannelCount': None,
    'EMGChannelCount': None,
    'MiscChannelCount': None,
    'TriggerChannelCount': None,
    'PowerLineFrequency': None,
    'EEGPlacementScheme': "10 percent system [need to confirm]",
    'EEGReference': "Common Average Reference",
    'EEGGround': "placed on Afz",
    'SoftwareFilters': "[Need to ask Tony]",
    'HardwareFilters': "[Need to ask Tony]",
    'RecordingDuration': None,
    'RecordingType': None
}



### Get list of files to convert ###

all_raw = []

for f in os.listdir(raw_data_path):
    
    # Ensure file is a valid EEG data file with a valid name
    study_id = re.findall(id_num_pattern, f)
    if not (raw_data_ext in f and len(study_id)):
        continue
        
    # Add padded study id and full path to file to list
    study_id = study_id[0].zfill(3)
    filepath = os.path.join(raw_data_path, f)
    all_raw.append((study_id, filepath))
    
    
### Initialize BIDS directory and process files ###

if not os.path.exists(bids_root):
    os.mkdir(bids_root)

print("\n###########################################")    
print("### Converting {0} files to BIDS format ###".format(len(all_raw)))
print("###########################################")
    
for study_id, filepath in all_raw:
    
    # Initialize path and read in data to MNE
    print("\n\n=== Converting '{0}'... ===\n".format(filepath))
    bids_path = BIDSPath(subject=study_id, task=taskname, root=bids_root)
    dat = mne.io.read_raw_edf(filepath, preload=False)
    
    # Fix channel types and round numbers before export
    dat.set_channel_types(channel_type_overrides)
    dat.info['sfreq'] = round(dat.info['sfreq'], 8)
    dat.info['lowpass'] = round(dat.info['lowpass'], 8)
    dat.info['highpass'] = round(dat.info['highpass'], 8)
    dat.info['line_freq'] = 60
    
    # Extract trigger event data from EEG annotations
    try:
        annot = mne.read_annotations(filepath)
        if len(annot) > 10:
            dat.set_annotations(annot)
            events, e_id = mne.events_from_annotations(dat, event_id=event_map)
            orig_time = dat.annotations.orig_time
        else:
            events = mne.find_events(dat, shortest_event=1, mask=65280,
                                     mask_type="not_and")
            orig_time = dat.info['meas_date']
        events = mne.pick_events(events, include=list(event_map.values()))
        annot_new = mne.annotations_from_events(
            events=events, sfreq=dat.info['sfreq'], orig_time=orig_time,
            event_desc=event_name_map, verbose=False
        )
        dat.set_annotations(annot_new)
    except (ValueError, RuntimeError):
        print("   * Unable to find any valid triggers, skipping...\n")
        continue
            
    # Acutally write out BIDS data
    write_raw_bids(dat, bids_path, verbose=False)
                   
    # Update sidecar files with correct metadata
    json_path = BIDSPath(subject=study_id, task=taskname, suffix='eeg',
                         extension='.json', root=bids_root)
    with open(json_path.fpath, 'r') as tmp_f:
        sidecar_json = json.load(tmp_f)
    file_info = metadata.copy()
    file_info.update(**sidecar_json)
    for k in metadata.keys():
        if metadata[k] != None:
            file_info[k] = metadata[k]
    _write_json(json_path.fpath, file_info, overwrite=True, verbose=False)
    

print("\n\n#################################")    
print("### BIDS conversion complete! ###")
print("#################################\n")                   
