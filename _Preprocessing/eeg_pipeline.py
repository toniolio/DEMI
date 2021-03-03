#########################################
## TraceLab EEG Preprocessing Pipeline ##
#########################################

# By: Austin Hurst

# Based on the 'run_full_prep.py' example from PyPREP (MIT Licence)


### Import required libraries #################################################

import os
import csv
import shutil
import random
import binascii
from collections import Counter

import mne
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from mne.preprocessing import ICA
from mne_bids import BIDSPath, read_raw_bids
from pyprep.prep_pipeline import PrepPipeline
from pyprep.find_noisy_channels import NoisyChannels
from pyprep.removeTrend import removeTrend
from pyprep.reference import Reference

from utils.save_edf import write_mne_edf


### Set general parameters ####################################################

matplotlib.use('Qt5Agg') # Only needed to avoid issues with ptpython
mne.set_log_level("WARNING")
seed = 530453080 # arbitrary, needed for results to be consistent between runs

task = 'tracelab'
datatype = 'eeg'
bids_root = 'bids'
output_root = 'output'

outdir = os.path.join(output_root, 'eeg')
plotdir = os.path.join(output_root, 'plots')
badsdir = os.path.join(outdir, 'bad')
noisy_bad_dir = os.path.join(badsdir, 'too_noisy')
ica_err_dir = os.path.join(badsdir, 'ica_err')
info_file = os.path.join(output_root, 'prep_info.csv')

perform_csd = True
interpolate_bads = True
max_interpolated = 0.25 # 25%
outfile_fmt = 'sub-{0}_eeg_prepped.edf'


def save_psd_plot(id_num, suffix, path, dat):

    plot_path = os.path.join(path, "sub-{0}_{1}.pdf".format(id_num, suffix))
    psd = dat.plot_psd(fmax = 100, show = False)
    psd.set_size_inches(6, 3)
    psd.savefig(plot_path, bbox_inches = 'tight')
    plt.close()


def save_channel_plot(id_num, suffix, path, dat):

    plot_path = os.path.join(path, "sub-{0}_{1}.png".format(id_num, suffix))
    ch_plot = dat.plot(
		n_channels=37, 
		duration=30, 
		scalings={"eeg": 1.5e-4, "eog": 5e2, "emg": 25e2},
		start=500
	)
    ch_plot.savefig(plot_path, bbox_inches = 'tight')
    plt.close()


def save_ica_plots(id_num, path, dat, ica, eog_scores):

    plot_path = os.path.join(path, "sub-{0}_ica_scores.pdf".format(id_num))
    ica_plot1 = ica.plot_scores(eog_scores, show = False)
    ica_plot1.savefig(plot_path, bbox_inches = 'tight')
    plt.close()

    ica_plot2 = ica.plot_properties(dat, picks = ica.exclude, show = False)
    n = 1
    for p in ica_plot2:
        plot_path = os.path.join(path, "sub-{0}_ica_{1}.pdf".format(id_num, n))
        p.savefig(plot_path, bbox_inches = 'tight')
        plt.close()
        n += 1

    plot_path = os.path.join(path, "sub-{0}_ica_sources.png".format(id_num))
    ica_plot3 = ica.plot_sources(dat, show = False)
    ica_plot3.savefig(plot_path, bbox_inches = 'tight')
    plt.close()


def save_bad_fif(dat, id_num, outdir_bad):
    if not os.path.isdir(outdir_bad):
        os.makedirs(outdir_bad)
    outpath = os.path.join(outdir_bad, outfile_fmt.format(id_num))
    outpath = outpath.replace('.edf', '.fif')
    dat.save(outpath)


def preprocess_eeg(id_num, random_seed=None):
    
    # Set important variables
    bids_path = BIDSPath(id_num, task=task, datatype=datatype, root=bids_root)
    plot_path = os.path.join(plotdir, "sub_{0}".format(id_num))
    if os.path.exists(plot_path):
        shutil.rmtree(plot_path)
    os.mkdir(plot_path)
    if not random_seed:
        random_seed = int(binascii.b2a_hex(os.urandom(4)), 16)
    random.seed(random_seed)
    id_info = {"id": id_num, "random_seed": random_seed}
    

    ### Load and prepare EEG data #############################################
    
    header = "### Processing sub-{0} (seed: {1}) ###".format(id_num, random_seed)
    print("\n" + "#" * len(header))
    print(header)
    print("#" * len(header) + "\n")

    # Load EEG data
    raw = read_raw_bids(bids_path, verbose=True)

    # Check if recording is complete
    complete = len(raw.annotations) >= 600

    # Add a montage to the data
    montage_kind = "standard_1005"
    montage = mne.channels.make_standard_montage(montage_kind)
    mne.datasets.eegbci.standardize(raw)
    raw.set_montage(montage)

    # Extract some info
    eeg_index = mne.pick_types(raw.info, eeg=True, eog=False, meg=False)
    ch_names = raw.info["ch_names"]
    ch_names_eeg = list(np.asarray(ch_names)[eeg_index])
    sample_rate = raw.info["sfreq"]

    # Make a copy of the data
    raw_copy = raw.copy()
    raw_copy.load_data()

    # Trim duplicated data (only needed for sub-005)
    annot = raw_copy.annotations
    file_starts = [a for a in annot if a['description'] == "file start"]
    if len(file_starts):
        duplicate_start = file_starts[0]['onset']
        raw_copy.crop(tmax=duplicate_start)

    # Make backup of EOG and EMG channels to re-append after PREP
    raw_other = raw_copy.copy()
    raw_other.pick_types(eog=True, emg=True, stim=False)

    # Prepare copy of raw data for PREP
    raw_copy.pick_types(eeg=True)

    # Plot data prior to any processing
    if complete:
        save_psd_plot(id_num, "psd_0_raw", plot_path, raw_copy)
        save_channel_plot(id_num, "ch_0_raw", plot_path, raw_copy)


    ### Clean up events #######################################################
    
    print("\n\n=== Processing Event Annotations... ===\n")

    event_names = [
        "stim_on", "red_on", "trace_start", "trace_end",
        "accuracy_submit", "vividness_submit"
    ]
    doubled = []
    wrong_label = []
    new_onsets = []
    new_durations = []
    new_descriptions = []

    # Find and flag any duplicate triggers
    annot = raw_copy.annotations
    trigger_count = len(annot)
    for i in range(1, trigger_count-1):
        a = annot[i]
        on_last = i+1 == trigger_count
        prev_trigger = annot[i-1]['description']
        next_onset = annot[i+1]['onset'] if not on_last else a['onset'] + 100
        # Determine whether duplicates are doubles or mislabeled
        if a['description'] == prev_trigger:
            if (next_onset - a['onset']) < 0.002:
                doubled.append(a)
            else:
                wrong_label.append(a)

    # Rename annotations to have meaningful names & fix duplicates
    for a in raw_copy.annotations:
        if a in doubled or a['description'] not in event_names:
            continue
        if a in wrong_label:
            index = event_names.index(a['description'])
            a['description'] = event_names[index + 1]
        new_onsets.append(a['onset'])
        new_durations.append(a['duration'])
        new_descriptions.append(a['description'])

    # Replace old annotations with new fixed ones
    if len(annot):
        new_annot = mne.Annotations(
            new_onsets, new_durations, new_descriptions,
            orig_time=raw_copy.annotations[0]['orig_time']
        )
        raw_copy.set_annotations(new_annot)

    # Check annotations to verify we have equal numbers of each
    orig_counts = Counter(annot.description)
    counts = Counter(raw_copy.annotations.description)
    print("Updated Annotation Counts:")
    for a in event_names:
        out = " - '{0}': {1} -> {2}"
        print(out.format(a, orig_counts[a], counts[a]))
    
    # Get info
    id_info['annot_doubled'] = len(doubled)
    id_info['annot_wrong'] = len(wrong_label)
    
    count_vals = [n for n in counts.values() if n != counts['vividness_submit']]
    id_info['equal_triggers'] = all(x == count_vals[0] for x in count_vals)
    id_info['stim_on'] = counts['stim_on']
    id_info['red_on'] = counts['red_on']
    id_info['trace_start'] = counts['trace_start']
    id_info['trace_end'] = counts['trace_end']
    id_info['acc_submit'] = counts['accuracy_submit']
    id_info['vivid_submit'] = counts['vividness_submit']

    if not complete:
        remaining_info = {
            'initial_bad': "NA", 'num_initial_bad': "NA",
            'interpolated': "NA", 'num_interpolated': "NA",
            'remaining_bad': "NA", 'num_remaining_bad': "NA"
        }
        id_info.update(remaining_info)
        e = "\n\n### Incomplete recording for sub-{0}, skipping... ###\n\n"
        print(e.format(id_num))
        return id_info


    ### Run components of PREP manually #######################################
    
    print("\n\n=== Performing CleanLine... ===")
    
    # Try to remove line noise using CleanLine approach
    linenoise = np.arange(60, sample_rate / 2, 60)
    EEG_raw = raw_copy.get_data() * 1e6
    EEG_new = removeTrend(EEG_raw, sample_rate=raw.info["sfreq"])
    EEG_clean = mne.filter.notch_filter(
        EEG_new,
        Fs=raw.info["sfreq"],
        freqs=linenoise,
        filter_length="10s",
        method="spectrum_fit",
        mt_bandwidth=2,
        p_value=0.01,
    )
    EEG_final = EEG_raw - EEG_new + EEG_clean
    raw_copy._data = EEG_final * 1e-6
    del linenoise, EEG_raw, EEG_new, EEG_clean, EEG_final
    
    # Plot data following cleanline
    save_psd_plot(id_num, "psd_1_cleanline", plot_path, raw_copy)
    save_channel_plot(id_num, "ch_1_cleanline", plot_path, raw_copy)

    # Perform robust re-referencing
    prep_params = {
        "ref_chs": ch_names_eeg,
        "reref_chs": ch_names_eeg
    }
    reference = Reference(
        raw_copy,
        prep_params,
        ransac=True,
        random_state=random_seed
    )
    print("\n\n=== Performing Robust Re-referencing... ===\n")
    reference.perform_reference()

    # If not interpolating bad channels, use pre-interpolation channel data
    if not interpolate_bads:
        reference.raw._data = reference.EEG_before_interpolation * 1e-6
        reference.interpolated_channels = []
        reference.still_noisy_channels = reference.bad_before_interpolation
        reference.raw.info["bads"] = reference.bad_before_interpolation

    # Plot data following robust re-reference
    save_psd_plot(id_num, "psd_2_reref", plot_path, reference.raw)
    save_channel_plot(id_num, "ch_2_reref", plot_path, reference.raw)

    # Re-append removed EMG/EOG/trigger channels
    raw_prepped = reference.raw.add_channels([raw_other])
    
    # Get info
    initial_bad = reference.noisy_channels_original["bad_all"]
    id_info['initial_bad'] = " ".join(initial_bad)
    id_info['num_initial_bad'] = len(initial_bad)
    
    interpolated = reference.interpolated_channels
    id_info['interpolated'] = " ".join(interpolated)
    id_info['num_interpolated'] = len(interpolated)
    
    remaining_bad = reference.still_noisy_channels
    id_info['remaining_bad'] = " ".join(remaining_bad)
    id_info['num_remaining_bad'] = len(remaining_bad)

    # Print re-referencing info
    print("\nRe-Referencing Info:")
    print(" - Bad channels original: {0}".format(initial_bad))
    if interpolate_bads:
        print(" - Bad channels after re-referencing: {0}".format(interpolated))
        print(" - Bad channels after interpolation: {0}".format(remaining_bad))
    else:
        print(" - Bad channels after re-referencing: {0}".format(remaining_bad))

    # Check if too many channels were interpolated for the participant
    prop_interpolated = len(reference.interpolated_channels) / len(ch_names_eeg)
    e = "### NOTE: Too many interpolated channels for sub-{0} ({1}) ###"
    if max_interpolated < prop_interpolated:
        print("\n")
        print(e.format(id_num, len(reference.interpolated_channels)))
        print("\n")


    ### Filter data and apply ICA to remove blinks ############################

    # Apply highpass & lowpass filters
    print("\n\n=== Applying Highpass & Lowpass Filters... ===")
    raw_prepped.filter(1.0, 50.0, fir_design='firwin')

    # Plot data following frequency filters
    save_psd_plot(id_num, "psd_3_filtered", plot_path, raw_prepped)
    save_channel_plot(id_num, "ch_3_filtered", plot_path, raw_prepped)

    # Perform ICA using EOG data on eye blinks
    print("\n\n=== Removing Blinks Using ICA... ===\n")
    ica = ICA(n_components=20, random_state=random_seed, method='picard')
    ica.fit(raw_prepped, decim=5)
    eog_indices, eog_scores = ica.find_bads_eog(raw_prepped)
    ica.exclude = eog_indices

    if not len(ica.exclude):
        err = " - Encountered an ICA error for sub-{0}, skipping for now..."
        print("\n")
        print(err.format(id_num))
        print("\n")
        save_bad_fif(raw_prepped, id_num, ica_err_dir)
        return id_info

    # Plot ICA info & diagnostics before removing from signal
    save_ica_plots(id_num, plot_path, raw_prepped, ica, eog_scores)

    # Remove eye blink independent components based on ICA
    ica.apply(raw_prepped)

    # Plot data following ICA
    save_psd_plot(id_num, "psd_4_ica", plot_path, raw_prepped)
    save_channel_plot(id_num, "ch_4_ica", plot_path, raw_prepped)


    ### Compute Current Source Density (CSD) estimates ########################

    if perform_csd:
        print("\n")
        print("=== Computing Current Source Density (CSD) Estimates... ===\n")
        raw_prepped = mne.preprocessing.compute_current_source_density(
            raw_prepped.drop_channels(remaining_bad)
        )

        # Plot data following CSD
        save_psd_plot(id_num, "psd_5_csd", plot_path, raw_prepped)
        save_channel_plot(id_num, "ch_5_csd", plot_path, raw_prepped)


    ### Write preprocessed data to new EDF ####################################

    if max_interpolated < prop_interpolated:
        if not os.path.isdir(noisy_bad_dir):
            os.makedirs(noisy_bad_dir)
        outpath = os.path.join(noisy_bad_dir, outfile_fmt.format(id_num))
    else:
        outpath = os.path.join(outdir, outfile_fmt.format(id_num))
    write_mne_edf(outpath, raw_prepped)
    
    print("\n\n### sub-{0} complete! ###\n\n".format(id_num))
    
    return id_info
    


### Run loop for all IDs ######################################################

if not os.path.isdir(output_root):
    os.mkdir(output_root)

if not os.path.isdir(outdir):
    os.mkdir(outdir)

if not os.path.isdir(plotdir):
    os.mkdir(plotdir)

ids = []
for sub in os.listdir(bids_root):
    if "sub-" in sub:
        ids.append(sub.split("-")[1])

processed_ids = []
if not os.path.isfile(info_file):
    open(info_file, 'w').close()
elif os.path.getsize(info_file) > 0:
    with open(info_file, 'r') as f:
        info_csv = csv.DictReader(f)
        processed_ids = [row['id'] for row in info_csv]

first_id = len(processed_ids) == 0
with open(info_file, 'a', newline='') as outfile:
    writer = csv.writer(outfile)
    print("")
    for sub in ids:
        if sub in processed_ids:
            print(" - sub-{0} already processed, skipping...\n".format(sub))
            continue
        info = preprocess_eeg(sub, random_seed=seed)
        if first_id:
            writer.writerow(list(info.keys()))
            first_id = False
        writer.writerow(list(info.values()))
