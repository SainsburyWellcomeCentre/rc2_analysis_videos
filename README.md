# RC2 Analysis of Videos

This project provides a suite of MATLAB scripts for analyzing facial and body motion energy, saccades, and pupil diameter for the translocator project's video data. The analyses target the effects of these behavioral metrics on neural firing rates, using extracted features from facial and body videos recorded during experimental trials.

This codebase is provided as a reference for the publication Velez-Fort 2024 and has been thoroughly documented. If you have any questions, please open an issue.

## Main Analysis Scripts

For the final analyses, the primary scripts are:

1. **`analyze_facial_motion_energy_and_firing_rate.m`**  
   Analyzes neuronal firing rate changes in relation to **facial motion energy (ME)** across different conditions. This script applies a facial ME mask and compares conditions with and without facial motion, providing key visualizations and summary statistics.

2. **`body_motion_energy_and_saccade_analysis.m`**  
   Examines **body motion energy** and **saccade events**, analyzing their impact on firing rates across conditions. Unity plots visualize differences in firing rates with and without body ME filtering, while modulation indices quantify population-level responses.

3. **`analyze_pupil_data.m`**  
   Processes **pupil diameter data** to evaluate its relationship with neural firing rates, calculating metrics such as modulation index.

## Supporting Data Acquisition Scripts

The following scripts assist in extracting essential data for the main analyses:

- **`extract_pupil_diameter.m`**: Extracts pupil diameter measurements from DeepLabCut outputs for downstream analysis.
- **`pupil_extract_saccade_frames.m`**: Identifies saccade events in video data by extracting frames from DeepLabCut outputs where saccades occur.

## Motion Energy Computation

Body and facial motion energy were precomputed using the [`rc2_analysis`](https://github.com/SainsburyWellcomeCentre/rc2_analysis) (see [`CameraProcessingHelper`](https://github.com/SainsburyWellcomeCentre/rc2_analysis/blob/8bea4a318d5af123b6086fa37a81fd2dcc9064ef/lib/classes/preprocess/CameraProcessingHelper.m)).

## Data Exploration Scripts

Additional exploratory scripts are included for investigating distributions, visualizing data, and testing various thresholds and parameters. These scripts provide preliminary insights and support hypothesis development but are not part of the final analysis pipeline.
