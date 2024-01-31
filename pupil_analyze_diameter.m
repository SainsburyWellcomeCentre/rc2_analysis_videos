% extracts frame numbers on which saccades occur from DeepLabCut tracks
%
%   if `save_csv` is true then a csv will be saved containing the frame
%   numbers of the saccades and the direction of the saccade.

session_id          = 'CAA-1112224_rec1';
fname               = ctl.file.camera0_dlc_pupil_slow(session_id);

% session_id          = 'CAA-1112874_rec1_001';
% fname               = ctl.file.camera0_dlc_pupil_slow(session_id);

lh_threshold        = 0.7;  % DLC likelihood threshold
fs                  = 60;   % sampling rate of camera
n_sds               = 4;    % number of S.D.s to use for detecting events
min_event_interval  = 0.2;  % seconds, remove events which occur within this time
save_csv            = true; % whether to save the saccades .csv file


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

% delta x and delta y
pupil = ([nasal.x, nasal.y] - [temporal.x, temporal.y]);

% remove uncertain values
valid = nasal.lh > lh_threshold & temporal.lh > lh_threshold;
pupil(~valid, :) = nan;

% % timebase of the traces
timebase        = (0:size(pupil, 1)-1) * dt;
frame_n         = 1:size(pupil, 1);

% diameter of the pupil
diameter = sqrt(pupil(:, 2).^2 + pupil(:, 1).^2);

%% plot pupil position and events in different colours

diff_centre = 392.5;

figure;
h_ax = axes();
hold on;

plot(h_ax, timebase, pupil(:, 1));
plot(h_ax, timebase, pupil(:, 2));
plot(h_ax, timebase, diameter);


xlabel(h_ax, 'Time (s)');
ylabel(h_ax, 'Pupil diameter (pixels)');

if save_csv
    writematrix(diameter, ctl.file.pupil_diameter_slow(session_id));
end
    



          
          
