% extract_saccade_frames.m
%
% Description:
%   This script extracts frame numbers on which saccades (rapid eye movements) occur 
%   from DeepLabCut (DLC) tracked pupil data. It calculates the pupil position and 
%   diameter over time and, if `save_csv` is set to true, saves the results in a .csv file.
%
% Inputs:
%   - session_id (string): Unique identifier for the recording session.
%   - fname (string): File path to the DLC-tracked pupil data for the session.
%   - lh_threshold (numeric): Threshold for DLC likelihood to filter out low-confidence points.
%   - fs (numeric): Sampling rate of the camera in Hz.
%   - n_sds (numeric): Standard deviation multiplier for detecting significant events.
%   - min_event_interval (numeric): Minimum interval in seconds between consecutive events.
%   - save_csv (boolean): Flag indicating whether to save the saccades data in a .csv file.
%
% Outputs:
%   - (If save_csv = true) .csv file containing frame numbers and directions of detected saccades.
%
% Usage:
%   Run this script after setting up session_id and other parameters.
%   Check and adjust `lh_threshold` and `n_sds` as necessary based on data quality.

session_id          = 'CAA-1112224_rec1';
fname               = ctl.file.camera0_dlc_pupil_slow(session_id); % Path to DLC tracking file

% Parameters for session data
lh_threshold        = 0.7;  % Likelihood threshold to filter out low-confidence data points
fs                  = 60;   % Sampling rate of the camera (frames per second)
n_sds               = 4;    % Multiplier for standard deviation to detect significant events
min_event_interval  = 0.2;  % Minimum interval between events, in seconds
save_csv            = true; % Whether to save saccade data to .csv

%% Read the DLC .csv and rename variables for clarity
tbl         = readtable(fname);  % Read DLC-tracked data file

% Extract and rename nasal and temporal tracking data (x, y, likelihood)
nasal.x     = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000;
nasal.y     = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_1;
nasal.lh    = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_2;

temporal.x  = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_3;
temporal.y  = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_4;
temporal.lh = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_5;

%% Calculate pupil position and diameter
dt = 1/fs;  % Time increment per frame

% Compute relative position between nasal and temporal points
pupil = ([nasal.x, nasal.y] - [temporal.x, temporal.y]);

% Filter out frames with low likelihood values (low confidence)
valid = nasal.lh > lh_threshold & temporal.lh > lh_threshold;
pupil(~valid, :) = nan;

% Define timebase for plotting and frame indexing
timebase        = (0:size(pupil, 1)-1) * dt;  % Time array for each frame
frame_n         = 1:size(pupil, 1);           % Frame numbers

% Calculate pupil diameter as Euclidean distance between nasal and temporal points
diameter = sqrt(pupil(:, 2).^2 + pupil(:, 1).^2);

%% Plot pupil position and diameter
diff_centre = 392.5;  % Centre offset for visualization (optional adjustment)

figure;
h_ax = axes();
hold on;

% Plot x, y components and diameter of pupil over time
plot(h_ax, timebase, pupil(:, 1), 'DisplayName', 'Pupil X Position');
plot(h_ax, timebase, pupil(:, 2), 'DisplayName', 'Pupil Y Position');
plot(h_ax, timebase, diameter, 'DisplayName', 'Pupil Diameter');

xlabel(h_ax, 'Time (s)');
ylabel(h_ax, 'Pupil Measurements (pixels)');
legend(h_ax, 'show');  % Display legend for each plotted line

%% Save data if save_csv is enabled
if save_csv
    % Save calculated pupil diameter to a CSV file
    writematrix(diameter, ctl.file.pupil_diameter_slow(session_id));
end
