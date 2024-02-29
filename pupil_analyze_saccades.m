% analyze firing rate around saccades

session_ids = {'CAA-1110262_rec1_001', ...
               'CAA-1110264_rec1_001', ...
               'CAA-1110265_rec1_001', ...
               'CAA-1112224_rec1_001'};

ctl = RC2Analysis();

camera_fs           = 60;
frames_around_event = (-30:30) ;
baseline_frames     = -5:1;

for ii = 1 : length(session_ids)
    
    % get probe ID from session ID
    probe_id = ctl.get_probe_id_from_session_id(session_ids{ii});
    
    % get formatted data
    data = ctl.load_formatted_data(probe_id);
    
    % get RVT trials
    trials = data.get_trials_with_trial_group_label('RVT');
    
    % get table of saccades
    saccade_tbl = ctl.load.camera0_saccades(session_ids{ii});
    
    % get the session object for the session ID
    session = data.get_session_with_id(session_ids{ii});
    
    % get the VISp cluster object
    clusters = data.VISp_clusters();
    
    % get pupil position from DeepLabCut .csv
    dlc_fname = ctl.file.camera0_dlc_pupil_slow(session_ids{ii});
    tbl = readtable(dlc_fname);
    nasal.x     = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000;
    temporal.x  = tbl.DLC_resnet50_pupil_trackingJan12shuffle1_350000_3;
    pupil = (nasal.x + temporal.x) / 2;
    
    % for temporal and nasal directions
    for kk = 1 : 2
        
        % get the frames
        if kk == 1
            str = 'nasal';
            idx = strcmp(saccade_tbl.direction, str);
            col = 'b';
            sd_col = [0.2, 0.2, 0.8];
        else
            str = 'temporal';
            idx = strcmp(saccade_tbl.direction, str);
            col = 'r';
            sd_col = [0.8, 0.2, 0.2];
        end
        
        % index of nasal or temporal saccades
        saccade_frames = saccade_tbl.saccade_frame(idx);
        
        fr = [];
        for clust_i = 1 : length(clusters)
            
            raster = RasterData(clusters(clust_i));
            raster.limits = [-0.5, 0.5];
            raster.trigger_times = session.camera_t(saccade_frames);
            
            % get firing rate average
            fr(:, clust_i) = raster.spike_convolution_avg();
        end
        
        fr_avg = mean(fr, 2);
        fr_sem = std(fr, [], 2) ./ sqrt(sum(~isnan(fr), 2));
        
        % get saccade events
        events = [];
        for frame_i = 1 : length(saccade_frames)
            idx = saccade_frames(frame_i) + frames_around_event;
            idx(idx > size(pupil, 1)) = size(pupil, 1);
            baseline_idx = saccade_frames(frame_i) + baseline_frames;
            events(:, frame_i) = pupil(idx, 1) - mean(pupil(baseline_idx, 1));
        end
        
        saccade_avg = mean(events, 2);
        saccade_sem = std(events, [], 2) ./ sqrt(sum(~isnan(events), 2));
        
        % find the trial (if it exists) for each trigger time
        associated_trial = nan(1, length(saccade_frames));
        
        for frame_i = 1 : length(saccade_frames)
            
            time_to_search = raster.trigger_times(frame_i);
            
            for trial_i = 1 : length(trials)
                if time_to_search > trials{trial_i}.probe_t(1) && time_to_search < trials{trial_i}.probe_t(end)
                    associated_trial(frame_i) = trial_i;
                    continue
                end
            end
        end
        
        saccades_in_trials = ~isnan(associated_trial);
        
        these_trials = trials(associated_trial(saccades_in_trials));
        these_times = raster.trigger_times(saccades_in_trials);
%         [motion_traces, motion_timebase] = data.get_traces_around_times(these_trials, these_times, [-0.2, 0.2], 10e3);
        
%         motion_avg = mean(motion_traces, 2);
%         motion_sem = std(motion_traces, [], 2) ./ sqrt(sum(~isnan(motion_traces), 2));
        
        figure;
        
        % plot saccade average
        h_ax = subplot(3, 1, 1);
        hold on;
        plot(h_ax, frames_around_event * (1/camera_fs), saccade_avg, 'color', col)
        plot(h_ax, frames_around_event * (1/camera_fs), saccade_avg + saccade_sem, 'color', sd_col)
        plot(h_ax, frames_around_event * (1/camera_fs), saccade_avg - saccade_sem, 'color', sd_col)
        line(h_ax, [0, 0], get(h_ax, 'ylim'), 'color', 'k');
        ylabel('Pixels');
        title(sprintf('%s saccade', str));
        
        % plot firing rate average
        h_ax = subplot(3, 1, 2);
        
        hold on;
        plot(h_ax, raster.common_t, fr_avg, 'color', col)
        plot(h_ax, raster.common_t, fr_avg+fr_sem, 'color', sd_col)
        plot(h_ax, raster.common_t, fr_avg-fr_sem, 'color', sd_col)
        
        h_ax.YLim = [0, max(h_ax.YLim)];
        line(h_ax, [0, 0], get(h_ax, 'ylim'), 'color', 'k');
        ylabel('Hz');
        title('Firing rate');
        
        % plot motion
%         h_ax = subplot(3, 1, 3);
%         hold on;
%         plot(h_ax, motion_timebase, motion_avg, 'color', col)
%         plot(h_ax, motion_timebase, motion_avg + motion_sem, 'color', sd_col)
%         plot(h_ax, motion_timebase, motion_avg - motion_sem, 'color', sd_col)
%         line(h_ax, [0, 0], get(h_ax, 'ylim'), 'color', 'k');
%         title('Motion');
%         ylabel('cm/s');
%         xlabel('Time (s)');
        
        FigureTitle(gcf, sprintf('%s', session_ids{ii}));
    end
end
