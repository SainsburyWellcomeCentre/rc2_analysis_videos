% extracts frame numbers on which saccades occur from DeepLabCut tracks
%
%   if `save_csv` is true then a csv will be saved containing the frame
%   numbers of the saccades and the direction of the saccade.

session_id          = 'CAA-1110264_rec1_001';
fname               = ctl.file.camera0_dlc_pupil_slow(session_id);

lh_threshold        = 0.7;  % DLC likelihood threshold
fs                  = 60;   % sampling rate of camera
n_sds               = 4;    % number of S.D.s to use for detecting events
min_event_interval  = 0.2;  % seconds, remove events which occur within this time
save_csv            = false; % whether to save the saccades .csv file


%% read the DLC .csv and rename the variables

tbl         = readtable(fname);

nasal.x     = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000;
nasal.y     = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_1;
nasal.lh    = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_2;

temporal.x  = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_3;
temporal.y  = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_4;
temporal.lh = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_5;


%%

dt = 1/fs;

% position of middle of pupil
pupil = ([nasal.x, nasal.y] + [temporal.x, temporal.y]) / 2;

% remove uncertain values
valid = nasal.lh > lh_threshold & temporal.lh > lh_threshold;
pupil(~valid, :) = nan;

% timebase of the traces
timebase        = (0:size(pupil, 1)-1) * dt;
frame_n         = 1:size(pupil, 1);

% central difference
t_diff          = timebase(2:end-1);
frame_n_diff    = frame_n(2:end-1);
pupil_diff      = (pupil(3:end, 1) - pupil(1:end-2, 1))/(2*dt);

% take average and standard deviation of difference trace
mean_diff       = nanmean(pupil_diff);
sd_diff         = nanstd(pupil_diff);

% find where the difference trace goes negative (nasal) or positive
% (temporal)
nasal_fast_move     = pupil_diff < mean_diff - n_sds*sd_diff;
temporal_fast_move  = pupil_diff > mean_diff + n_sds*sd_diff;

%
nasal_fast_move_onset       = [false; diff(nasal_fast_move) == 1];
temporal_fast_move_onset    = [false; diff(temporal_fast_move) == 1];

% restrict to onset movements followed by another fast move
nasal_fast_move_onset       = [nasal_fast_move_onset(1:end-1) & nasal_fast_move(2:end); false];
temporal_fast_move_onset    = [temporal_fast_move_onset(1:end-1) & temporal_fast_move(2:end); false];

% get the frame on which the fast move occurs
nasal_fast_move_frames      = frame_n_diff(nasal_fast_move_onset);
temporal_fast_move_frames   = frame_n_diff(temporal_fast_move_onset);

% the times
nasal_times                 = timebase(nasal_fast_move_frames);
temporal_times              = timebase(temporal_fast_move_frames);


% remove events which happen within 200ms of each other
delta_times                 = nasal_times(:) - temporal_times(:)';
delta_times                 = abs(delta_times) < min_event_interval;

nasal_rm                    = sum(delta_times, 2) > 0;
temporal_rm                 = sum(delta_times, 1) > 0;

nasal_fast_move_frames(nasal_rm)        = [];
temporal_fast_move_frames(temporal_rm)  = [];
nasal_times(nasal_rm)                   = [];
temporal_times(temporal_rm)             = [];

fprintf('Removing %i/%i events within 200ms of each other\n', sum(nasal_rm) + sum(temporal_rm), length(nasal_rm)+length(temporal_rm));


%% plot pupil position and events in different colours

diff_centre = 392.5;

figure;
h_ax = axes();
hold on;

plot(h_ax, timebase, pupil(:, 1));

for ii = 1 : length(nasal_times)
    
    line(h_ax, nasal_times([ii, ii]), get(h_ax, 'ylim'), 'color', 'r');
end

for ii = 1 : length(temporal_times)
    
    line(h_ax, temporal_times([ii, ii]), get(h_ax, 'ylim'), 'color', 'g');
end

set(h_ax, 'plotboxaspectratio', [3, 1, 1])

plot(h_ax, t_diff, diff_centre + (pupil_diff - mean_diff) / (2*n_sds*sd_diff));

line(h_ax, get(h_ax, 'xlim'), diff_centre + [0.5, 0.5], 'color', 'k')
line(h_ax, get(h_ax, 'xlim'), diff_centre - [0.5, 0.5], 'color', 'k')

plot(h_ax, t_diff, diff_centre - 2.5 + nasal_fast_move);
plot(h_ax, t_diff, diff_centre - 1.5 + temporal_fast_move);

xlabel(h_ax, 'Time (s)');
ylabel(h_ax, 'Pupil position (pixels)');


%% plot pupil trace around events

frames_around_event = -30:30;
baseline_frames = -5:-1;

nasal_events = [];
for ii = 1 : length(nasal_fast_move_frames)
    
    idx = nasal_fast_move_frames(ii) + frames_around_event;
    baseline_idx = nasal_fast_move_frames(ii) + baseline_frames;
    nasal_events(:, ii) = pupil(idx, 1) - mean(pupil(baseline_idx, 1));
end

figure
h_ax = axes();
hold on;
plot(h_ax, frames_around_event, nasal_events);
line(h_ax, [0, 0], get(gca, 'ylim'));


temporal_events = [];
for ii = 1 : length(temporal_fast_move_frames)
    
    idx = temporal_fast_move_frames(ii) + frames_around_event;
    baseline_idx = temporal_fast_move_frames(ii) + baseline_frames;
    temporal_events(:, ii) = pupil(idx, 1) - mean(pupil(baseline_idx, 1));
end

figure
h_ax = axes();
hold on;
plot(h_ax, frames_around_event, temporal_events);
line(h_ax, [0, 0], get(gca, 'ylim'));
xlabel(h_ax, 'Frame #');
ylabel(h_ax, 'Pupil position (\Delta pixels)');

%% plot average pupil position around events

nasal_event_avg = nanmean(nasal_events, 2);
temporal_event_avg = nanmean(temporal_events, 2);

figure
h_ax = axes();
hold on;
plot(h_ax, frames_around_event, nasal_event_avg, 'color', 'b');
plot(h_ax, frames_around_event, temporal_event_avg, 'color', 'r');
line(h_ax, [0, 0], get(gca, 'ylim'));
xlabel(h_ax, 'Frame #');
ylabel(h_ax, 'Pupil position (\Delta pixels)');


%% write csv

n_nasal_events      = length(nasal_fast_move_frames);
n_temporal_events   = length(temporal_fast_move_frames);

all_frames          = [nasal_fast_move_frames(:); temporal_fast_move_frames(:)];
group               = [ones(n_nasal_events, 1); 2*ones(n_temporal_events, 1)];

[all_frames, sort_idx] = sort(all_frames);
group               = group(sort_idx);

% check
assert(isequal(all_frames(group == 1), nasal_fast_move_frames(:)), 'something wrong');
assert(isequal(all_frames(group == 2), temporal_fast_move_frames(:)), 'something wrong');

direction           = cell(length(group), 1);
for ii = 1 : length(group)
    if group(ii) == 1
        direction{ii} = 'nasal';
    elseif group(ii) == 2
        direction{ii} = 'temporal';
    end
end
        

% create table
tbl                 = table(all_frames, direction, 'VariableNames', {'saccade_frame', 'direction'});

% write to csv
if save_csv
    writetable(tbl, 'camera0_saccades.csv');
end
          
          
