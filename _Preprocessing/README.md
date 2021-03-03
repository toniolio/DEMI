# DEMI EEG Preprocessing Pipeline

## Processing Steps

1. Repair mislabeled or duplicated event triggers
2. Remove line noise from the data using CleanLine
3. Perform robust rereferencing and re-interpolate bad channels using PyPREP
4. Apply highpass (1 Hz) and lowpass (50 Hz) filters to the data
5. Use ICA and EOG channel data to detect and remove components of the signal associated with blinks
6. Save the preprocessed data files to EDF+

## Requirements

- Python 3.7 (newer versions may work, but are untested)
- A reasonably fast computer with at least 8 GB RAM
- The `pipenv` package

## Instructions

### Running the pipeline

1. First, copy or move the contents of the BIDS-formatted dataset to the project's `bids` folder (to convert the loose EDF data to an organized BIDS dataset, you can use the `edf2bids.py` script as described below).
2. Navigate to the project folder in a terminal window.
3. Run `pipenv install` at the terminal prompt (NOTE: if using a Python version other than 3.7, you may need to modify the project's Pipfile accordingly first). This will create a clean virtual environment for the project and install all the pipeline's dependencies within it.
4. When you're ready, run `pipenv run python eeg_pipeline.py` in the terminal to start the process.
5. Wait for several hours as the pipeline completes!

**NOTE:** Running the pipeline will consume much of your system memory and may make your computer slow or unresponsive, so it may be a good idea to run it overnight.

### Creating the BIDS dataset (optional)

This pipline includes a script, `edf2bids.py`, which was used to convert the loose EDF files for the project into the BIDS-formatted dataset used for the preprocessing script. This script does not modify the original EDF files at all, but instead reorganizes them and adds additional metadata files to make working with and sharing the dataset easier (e.g. properly labelling the channel types of the EMG and EOG channels so MNE doesn't interpret them as EEG, adding metadata about the amplifier and software used for data collection, etc.). This script is project-specific, but can likely be adapted easily for other projects.

To use the script, add all .edf files into a folder named "edfs" in the root of this environment, then run `pipenv run python edf2bids.py`. When complete, the `bids` folder should be populated with the resulting BIDS dataset.

At present, the `edf2bids.py` script ignores all split EDFs, as well as any EDFs that don't contain any triggers or annotations.

## Output

The pipeline produces three forms of output:

1. The preprocessed EEG data in EDF+ format, located in the `output/eeg` folder
2. Plots visualizing the effects of the preprocessing steps, located in the `output/plots` folder
3. Data on the preprocessing for all particpants, located in the `output/prep_info.csv` file

Sessions flagged as unusable due to excessive noise or problems with ICA removal of blinks are saved as MNE `.fif` files in the `output/bads/too_noisy` and `output/bads/ica_err` folders, respectively. This allows partially-processed problem files to be inspected manually if needed.
