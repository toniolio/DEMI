#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" pyedf is a python package to read from and write EEG data to European Data 
    Format files. Since EDF is such a widely used format, there exist multiple 
    Python implementations for reading and writing EDF files. However, most of 
    these Python modules consist of wrappers around the C-code implementation, 
    which makes installation more cumbersome and reduces portability. This 
    implementation is in pure python with limited dependencies on external 
    packages while having support for Python 2.7 and 3.x.

    Note: the EDF header is represented as a tuple of (meas_info, chan_info)
        meas_info should have the following: [
            'record_length', 'file_ver', 'hour', 'subject_id',
            'recording_id', 'n_records', 'month', 'subtype',
            'second', 'nchan', 'data_size', 'data_offset', 
            'lowpass', 'year', 'highpass', 'day', 'minute']
        chan_info should have the following: [
            'physical_min', 'transducers', 'physical_max', 
            'digital_max', 'ch_names', 'n_samps', 'units', 'digital_min']
    
    The EDF standard was introduced in 1992. The extension of EDF with 
    annotations was first described in 1998 and more formalized with the EDF+ 
    standard that was published in 2003. To learn more about both standards and 
    implementation details, check out https://www.edfplus.info/index.html
"""
__version__ = "0.2.0" # NOTE: Heavily modified by Austin Hurst
__author__ = "Robert Oostenveld <r.oostenveld@gmail.com> | Sandeepan B <bsandeepan95.work@gmail.com>"
__copyright__ = "https://bids.neuroimaging.io/"
__credits__ = ["various"]
__license__ = "BSD 3-Clause License"
__maintainer__ = "https://github.com/bids-standard"
__status__ = "Production"


import os
import datetime
from copy import deepcopy
from math import ceil, floor
from warnings import warn

import numpy as np


def is_numeric(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def padtrim(buf, size, pad = ' '):
    return str(buf).ljust(size, pad)[:size]

def validate_header_field(x, size, name):
    s = str(int(x)) if is_numeric(x) else str(x)
    if len(s) > size:
        typestr = "number" if is_numeric(x) else "string"
        if isinstance(x, (np.floating, float)):
            s = "{:.3f}".format(x).rstrip("0")
        e = ("The {0} '{1}' is too long to fit in the '{2}' field of "
             "the EDF+ header (maximum {3} characters).")
        raise ValueError(e.format(typestr, s, name, size))

def writebyte(file, content, encoding='ascii'):
    try:
        # Py3 onwards bytes and strings are separate data format
        content = bytes(content, encoding)
    except TypeError :
        # Py2.7 Support
        content = content.encode(encoding)
    except Exception as e:
        print(type(e))
        print(str(e))
        print(
            "If you see this message, please go to " + 
            "https://github.com/bids-standard/pyedf/issues" + 
            " and open an issue there regarding this. Thank you!")
    finally:
        file.write(content)

def set_offsets(chan_obj):
    offsets = []
    calibrations = []
    for i in range(len(chan_obj['ch_names'])):
        # Truncate floats for physical min/max to match header values
        pmin = float(padtrim(chan_obj['physical_min'][i], 8))
        pmax = float(padtrim(chan_obj['physical_max'][i], 8))
        dmin = chan_obj['digital_min'][i]
        dmax = chan_obj['digital_max'][i]
        # Calculate and append offsets and calibrations
        calibrate = (pmax - pmin) / (dmax - dmin)
        offset = (pmax / calibrate) - dmax
        offsets.append(offset)
        calibrations.append(calibrate)
    return (offsets, calibrations)

def verify_ascii(x):
    err = "'{0}' contains unsupported characters (only printable ASCII allowed)"
    vals = [ord(c) for c in x]
    if min(vals) == 32:
        raise ValueError("Header fields cannot contain spaces")
    elif min(vals) < 32 or max(vals) > 126:
        raise ValueError(err.format(x))

def date_to_str(dt, fmt="%d-%b-%Y"):
    try:
        datestr = dt.strftime(fmt)
    except AttributeError:
        raise ValueError("All dates must be passed as Python datetime objects")
    return datestr.upper()

def create_patient_id(pid="X", sex="X", birthdate="X", name="X"):
    # First, check if inputs are valid
    verify_ascii(pid)
    if sex.upper()[0] not in ["M", "F", "X"]:
        raise ValueError("EDF+ only supports 'M', 'F', and 'X' as valid sex values")
    sex = sex.upper()[0]
    if birthdate != "X":
        birthdate = date_to_str(birthdate)
    verify_ascii(name)
    # Once validated, create and return string
    outstr = " ".join([pid, sex, birthdate, name])
    if len(outstr) > 80:
        raise ValueError("Combined ID string cannot exceed 80 characters")
    return outstr

def create_recording_id(startdate="X", admin_code="X", technician="X", equipment="X"):
    # First, check if inputs are valid
    if startdate != "X":
        startdate = date_to_str(startdate)
    verify_ascii(admin_code)
    verify_ascii(technician)
    verify_ascii(equipment)
    # Once validated, create and return string
    outstr = " ".join(["Startdate", startdate, admin_code, technician, equipment])
    if len(outstr) > 80:
        raise ValueError("Combined ID string cannot exceed 80 characters")
    return outstr

def create_filter_str(highpass=None, lowpass=None, notch=None, more=None):
    outstr = ""
    if highpass:
        fmt = "{:.0f}" if highpass == int(highpass) else "{:.1f}"
        outstr += "HP:{0}Hz".format(fmt.format(highpass))
    if lowpass:
        pad = " " if len(outstr) else ""
        fmt = "{:.0f}" if lowpass == int(lowpass) else "{:.1f}"
        outstr += "{0}LP:{1}Hz".format(pad, fmt.format(lowpass))
    if notch:
        pad = " " if len(outstr) else ""
        fmt = "{:.0f}" if notch == int(notch) else "{:.1f}"
        outstr += "{0}N:{1}Hz".format(pad, fmt.format(notch))
    if more:
        verify_ascii(more)
        pad = " " if len(outstr) else ""
        outstr += (pad + more)
    return outstr

def create_annotation(onset, duration='', description=''):
    tmp = "+{0}\x15{1}\x14{2}\x14\x00"
    onset_fmt = "{:.0f}" if onset == int(onset) else "{:.3f}"
    onset_str = onset_fmt.format(onset)
    if duration != '':
        dur_fmt = "{:.0f}" if duration == int(duration) else "{:.3f}"
        duration = dur_fmt.format(duration)
    return tmp.format(onset_str, duration, description)



class EDFWriter():

    def __init__(self, fname=None):

        self._initial_state()
        if fname:
            self.open(fname)
    

    def _initial_state(self):

        self.fname = None
        self.meas_info = None
        self.chan_info = None
        self.calibrate = None
        self.offset    = None
        self.bdf = False
        self.n_records = 0
        self.annot_ch = 0
        self.annot_samps = 0
        self.annotations = []
    

    def open(self, fname):

        with open(fname, 'wb') as fid:
            assert fid.tell() == 0
        self.fname = fname


    def close(self):

        # Update the number of records in the header before closing
        with open(self.fname, 'r+b') as fid:
            assert fid.tell() == 0
            fid.seek(236)
            writebyte(fid, padtrim(self.n_records, 8)) 
            
        self._initial_state()


    def add_annotation(self, onset, duration, description):

        self.annotations.append([onset, duration, description])


    def write_header(self, header):

        meas_info = header[0]
        chan_info = header[1]

        with open(self.fname, 'wb') as fid:
            assert fid.tell() == 0

            # Initialize important info
            nchan = int(meas_info['nchan'])
            if meas_info['subtype'] in ('24BIT', 'bdf'):
                meas_info['data_size'] = 3  # 24-bit (3 byte) integers
                typestr = 'BDF+C'
                dmin, dmax = -8388608, 8388607
                self.bdf = True
            else:
                meas_info['data_size'] = 2  # 16-bit (2 byte) integers
                typestr = 'EDF+C'
                dmin, dmax = -32768, 32767

            # Fill in missing or incomplete information
            if not 'subject_id' in meas_info:
                meas_info['subject_id'] = create_patient_id()
            
            if not 'recording_id' in meas_info:
                meas_info['recording_id'] = create_recording_id()
            
            if not 'subtype' in meas_info:
                meas_info['subtype'] = 'edf'

            if not 'startdate' in meas_info:
                meas_info['startdate'] = '01.01.85'

            if not 'starttime' in meas_info:
                meas_info['starttime'] = '00.00.00'

            if not 'record_count' in meas_info:
                meas_info['record_count'] = -1
            
            if not 'ch_names' in chan_info:
                chan_info['ch_names'] = [str(i) for i in range(nchan)]

            self.annot_ch = 0
            annot_ch_info = {}
            annot_ch_name = typestr[0:3] + " Annotations"
            if not annot_ch_name in chan_info['ch_names']:
                self.annot_ch = 1
                self.annot_samps = max(100, chan_info['n_samps'][0])
                annot_ch_info = {
                    'ch_names': annot_ch_name, 'transducers': '', 'units': '',
                    'physical_min': 0, 'physical_max': 1,
                    'digital_min': dmin, 'digital_max': dmax,
                    'prefilters': '', 'n_samps': self.annot_samps
                }

            meas_size = 256
            chan_size = 256 * (nchan + self.annot_ch)

            # Write out non-channel metadata to the header
            writebyte(fid, (padtrim('0', 8)))
            writebyte(fid, (padtrim(meas_info['subject_id'], 80)))
            writebyte(fid, (padtrim(meas_info['recording_id'], 80)))

            writebyte(fid, (padtrim(meas_info['startdate'], 8)))
            writebyte(fid, (padtrim(meas_info['starttime'], 8)))
            
            writebyte(fid, (padtrim(meas_size + chan_size, 8)))
            writebyte(fid, (padtrim(typestr, 44)))

            # the final n_records should be inserted on byte 236
            writebyte(fid, (padtrim(meas_info['record_count'], 8)))  
            writebyte(fid, (padtrim(meas_info['record_length'], 8)))
            writebyte(fid, (padtrim(meas_info['nchan'] + self.annot_ch, 4)))

            # Write out per-channel metadata
            chan_fields = [
                'ch_names', 'transducers', 'units', 'physical_min', 'physical_max', 
                'digital_min', 'digital_max', 'prefilters', 'n_samps'
            ]
            bytes_per_field = {
                'ch_names': 16, 'transducers': 80, 'units': 8,
                'physical_min': 8, 'physical_max': 8,
                'digital_min': 8, 'digital_max': 8,
                'prefilters': 80, 'n_samps': 8
            }
            for key in chan_fields:
                # Double-check and sanitize field inputs
                if not key in chan_info:
                    if key in ['transducers', 'prefilters', 'units']:
                        chan_info[key] = [''] * nchan
                    else:
                        e = "Missing values for required channel info field '{0}'"
                        raise ValueError(e.format(key))
                if len(chan_info[key]) != nchan:
                    e1 = "Length of channel info field '{0}' does not match number of channels "
                    e2 = "(expected: {0}, recieved: {1})".format(nchan, len(chan_info[key]))
                    raise ValueError((e1.format(key) + e2))
                chan_info[key] = np.asarray(chan_info[key])
                # Actually write out data
                for val in chan_info[key]:
                    validate_header_field(val, bytes_per_field[key], key)
                    writebyte(fid, padtrim(val, bytes_per_field[key]))
                # If needed, add annotation channel info
                if self.annot_ch == 1:
                    writebyte(fid, padtrim(annot_ch_info[key], bytes_per_field[key]))

            for i in range(nchan + self.annot_ch):
                writebyte(fid, (' ' * 32)) # reserved
            
            # Set data offset
            meas_info['data_offset'] = fid.tell()

        self.meas_info = meas_info
        self.chan_info = chan_info
        self.offset, self.calibrate = set_offsets(chan_info)

        channels = list(range(nchan))
        for ch in channels:
            if self.calibrate[ch] < 0:
                self.calibrate[ch] = 1
                self.offset[ch] = 0



    def write_data(self, data):

        annot = deepcopy(self.annotations)
        nchan = self.meas_info['nchan']
        seconds_per_block = self.meas_info['record_length']

        # Verify equal numbers of blocks per channel
        num_blocks = []
        for i in range(nchan):
            num_blocks.append(int(len(data[i]) / self.chan_info['n_samps'][i]))
        if not all([n == num_blocks[0] for n in num_blocks]):
            raise ValueError("Uneven number of blocks per channel not supported.")

        with open(self.fname, 'ab') as fid:
            assert fid.tell() > 0
            for block in range(1, num_blocks[0] + 1):
                for i in range(nchan):

                    samps_per_block = self.chan_info['n_samps'][i]
                    start = (block - 1) * samps_per_block
                    end = block * samps_per_block
                    raw = deepcopy(data[i][start:end])

                    if min(raw) < self.chan_info['physical_min'][i]:
                        warn('Value exceeds physical_min: ' + str(min(raw)))
                    if max(raw) > self.chan_info['physical_max'][i]:
                        warn('Value exceeds physical_max: '+ str(max(raw)))

                    raw = (raw / self.calibrate[i]) - self.offset[i]

                    raw = np.asarray(raw).astype(np.int16)
                    raw.tofile(fid)

                if self.annot_ch == 1:
                    t = (block - 1) * seconds_per_block
                    t_next = block * seconds_per_block
                    block_annot = create_annotation(t)
                    while len(annot) and annot[0][0] <= t_next:
                        onset, dur, desc = annot.pop(0)
                        block_annot += create_annotation(onset, dur, desc)
                    writebyte(fid, (padtrim(block_annot, self.annot_samps * 2, pad='\x00')))

                self.n_records += 1
