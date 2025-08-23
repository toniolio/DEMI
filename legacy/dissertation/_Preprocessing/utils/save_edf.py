# -*- coding: utf-8 -*-

import os
import numpy as np
import utils.EDF as pyedf


def write_mne_edf(fname, raw):

    ext = os.path.splitext(fname)[-1]
    if ext == ".edf":
        filetype = 'edf'
        data_size = 2
        dmin, dmax = -32768, 32767
    elif ext == ".bdf":
        filetype = 'bdf'
        data_size = 3 # NOTE: BDF export isn't actually supported yet
        dmin, dmax = -8388608, 8388607
    
    data = raw.get_data()
    nchan = raw.info["nchan"]
    date = raw.info["meas_date"]
    default_lopass = (raw.info['sfreq'] / 2.0) - 1
    hipass = raw.info['highpass'] if raw.info['highpass'] > 0 else None
    lopass = raw.info['lowpass'] if raw.info['lowpass'] < default_lopass else None
    filt = pyedf.create_filter_str(hipass, lopass)

    ch_types = np.asarray(raw.get_channel_types())
    data[ch_types == 'eeg', :] *= 1e6  # scale normal EEG to microvolts
    data[ch_types == 'csd', :] *= 1e4  # scale CSD EEG to millivolts

    ch_units = []
    for ch_type in ch_types:
        if ch_type in ['eog', 'emg']:
            ch_units.append('V')
        elif ch_type == 'eeg':
            ch_units.append('uV')
        elif ch_type == 'csd':
            ch_units.append('mV/m^2')
        else:
            ch_units.append('')

    meas_info = {
        'subtype': filetype,
        'recording_id': pyedf.create_recording_id(startdate=date, technician="PyPREP"),
        'startdate': pyedf.date_to_str(date, fmt="%d.%m.%y"),
        'starttime': pyedf.date_to_str(date, fmt="%H.%M.%S"),
        'record_count': int(len(data[0]) / raw._raw_extras[0]['n_samps'][0]),
        'record_length': raw._raw_extras[0]['record_length'][0],
        'nchan': len(data)
    }
    chan_info = {
        'ch_names': raw.info["ch_names"],
        'transducers': [''] * nchan,
        'units': ch_units,
        'physical_min': data.min(axis=1),
        'physical_max': data.max(axis=1),
        'digital_min': [dmin] * nchan,
        'digital_max': [dmax] * nchan,
        'prefilters': [filt] * nchan,
        'n_samps':raw._raw_extras[0]['n_samps'][0:nchan]
    }

    f = pyedf.EDFWriter(fname)
    f.write_header([meas_info, chan_info])
    for ann in raw.annotations:
        f.add_annotation(ann["onset"], ann["duration"], ann["description"])
    f.write_data(data)
    f.close()
    