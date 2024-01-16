% Produce unity plots with no change, with removal of high motion energy
% periods (for repeat trials) and with removal of periods just after
% movements of the pupil.

experiment_groups           = 'visual_flow';
trial_types                 = {'RVT', 'RV', {'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};
duration_to_remove          = 0.25;     % duration after saccade to remove

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);

motion_fr_store             = cell(1, length(probe_ids));
stationary_fr_store         = cell(1, length(probe_ids));

motion_fr_store_pupil       = cell(1, length(probe_ids));
stationary_fr_store_pupil   = cell(1, length(probe_ids));

motion_fr_store_ME          = cell(1, length(probe_ids));
stationary_fr_store_ME      = cell(1, length(probe_ids));

stationary_time             = cell(1, length(probe_ids));
motion_time                 = cell(1, length(probe_ids));

stationary_time_pupil       = cell(1, length(probe_ids));
motion_time_pupil           = cell(1, length(probe_ids));

stationary_time_ME          = cell(1, length(probe_ids));
motion_time_ME              = cell(1, length(probe_ids));

for probe_i = 1 : length(probe_ids)
    
    data   = ctl.load_formatted_data(probe_ids{probe_i});
    clusters  = data.VISp_clusters();
    
    % get table of saccades
    sessions                        = data.motion_sessions;
    saccade_tbl                     = ctl.load.camera0_saccades(sessions{1}.session_id);
    saccade_times                   = sessions{1}.camera_t(saccade_tbl.saccade_frame);
    
    motion_fr_store{probe_i} = cell(1, length(trial_types));
    stationary_fr_store{probe_i} = cell(1, length(trial_types));
    
    motion_fr_store_pupil{probe_i} = cell(1, length(trial_types));
    stationary_fr_store_pupil{probe_i} = cell(1, length(trial_types));
    
    motion_fr_store_ME{probe_i} = cell(1, length(trial_types));
    stationary_fr_store_ME{probe_i} = cell(1, length(trial_types));
    
    stationary_time{probe_i} = cell(1, length(trial_types));
    motion_time{probe_i} = cell(1, length(trial_types));
    
    stationary_time_pupil{probe_i} = cell(1, length(trial_types));
    motion_time_pupil{probe_i} = cell(1, length(trial_types));
    
    stationary_time_ME{probe_i} = cell(1, length(trial_types));
    motion_time_ME{probe_i} = cell(1, length(trial_types));
    
    for type_i = 1 : length(trial_types)
        
        trials = data.get_trials_with_trial_group_label(trial_types{type_i});
        
        motion_fr_store{probe_i}{type_i} = nan(length(trials), length(clusters));
        stationary_fr_store{probe_i}{type_i} = nan(length(trials), length(clusters));
        
        motion_fr_store_pupil{probe_i}{type_i} = nan(length(trials), length(clusters));
        stationary_fr_store_pupil{probe_i}{type_i} = nan(length(trials), length(clusters));
        
        motion_fr_store_ME{probe_i}{type_i} = nan(length(trials), length(clusters));
        stationary_fr_store_ME{probe_i}{type_i} = nan(length(trials), length(clusters));
        
        stationary_time{probe_i}{type_i} = nan(length(trials), 1);
        motion_time{probe_i}{type_i} = nan(length(trials), length(clusters));
        
        stationary_time_pupil{probe_i}{type_i} = nan(length(trials), length(clusters));
        motion_time_pupil{probe_i}{type_i} = nan(length(trials), length(clusters));
        
        stationary_time_ME{probe_i}{type_i} = nan(length(trials), length(clusters));
        motion_time_ME{probe_i}{type_i} = nan(length(trials), length(clusters));
        
        for trial_i = 1 : length(trials)
            
            trial                = trials{trial_i}.to_aligned;
            
            % find the saccades in this trial
            timebase            = trial.probe_t;
            idx = saccade_times > timebase(1) & saccade_times < timebase(end);
            trial_saccade_times = saccade_times(idx);
            
            saccade_mask        = false(size(timebase));
            
            for sac_i = 1 : length(trial_saccade_times)
                saccade_mask = saccade_mask | ...
                                    (timebase > trial_saccade_times(sac_i) & timebase < ...
                                     trial_saccade_times(sac_i) + duration_to_remove);
            end
            
            
            if type_i > 2
                
                % get this trial and original trial
                original_trial              = trial.original_trial;
                
                % mask of motion in the original trial
                original_motion_mask        = original_trial.motion_mask;
                original_stationary_mask    = original_trial.stationary_mask;
                
                % camera motion in original trial
                original_motion_energy      = original_trial.camera1;
                
                % select the motion periods and get the lowest 5th prctile
                cam_motion_original         = original_motion_energy(original_motion_mask);
                motion_threshold            = prctile(cam_motion_original, 5);
                
                % camera motion in replay trial, aligned
                replay_motion_energy        = trial.camera1;
            end
            
            
            
            
            stationary_mask = trial.stationary_mask;
            motion_mask = trial.motion_mask;
            
            stationary_mask_pupil   = stationary_mask & ~saccade_mask;
            motion_mask_pupil       = motion_mask & ~saccade_mask;
            
            if type_i > 2
                stationary_mask_ME      = stationary_mask & (replay_motion_energy < motion_threshold);
                motion_mask_ME          = motion_mask & (replay_motion_energy < motion_threshold);
                stationary_time_ME{probe_i}{type_i}(trial_i) = sum(stationary_mask_ME)/10e3;
                motion_time_ME{probe_i}{type_i}(trial_i) = sum(motion_mask_ME)/10e3;
            end
            
            stationary_time{probe_i}{type_i}(trial_i) = sum(stationary_mask)/10e3;
            motion_time{probe_i}{type_i}(trial_i) = sum(motion_mask)/10e3;
            stationary_time_pupil{probe_i}{type_i}(trial_i) = sum(stationary_mask_pupil)/10e3;
            motion_time_pupil{probe_i}{type_i}(trial_i) = sum(motion_mask_pupil)/10e3;
            
            
            for clust_i = 1 : length(clusters)
                
                fr = clusters(clust_i).fr.get_convolution(trial.probe_t);
                
                motion_fr_store{probe_i}{type_i}(trial_i, clust_i) = mean(fr(motion_mask));
                stationary_fr_store{probe_i}{type_i}(trial_i, clust_i) = mean(fr(stationary_mask));
                
                motion_fr_store_pupil{probe_i}{type_i}(trial_i, clust_i) = mean(fr(motion_mask_pupil));
                stationary_fr_store_pupil{probe_i}{type_i}(trial_i, clust_i) = mean(fr(stationary_mask_pupil));
                
                if type_i > 2
                    motion_fr_store_ME{probe_i}{type_i}(trial_i, clust_i) = mean(fr(motion_mask_ME));
                    stationary_fr_store_ME{probe_i}{type_i}(trial_i, clust_i) = mean(fr(stationary_mask_ME));
                end
            end
        end
    end
end


%%
for type_i = 1 : length(trial_types)
    
    median_stationary_fr    = [];
    median_motion_fr        = [];
    p_val                   = [];
    direction               = [];
    
    median_stationary_fr_pupil    = [];
    median_motion_fr_pupil        = [];
    p_val_pupil                   = [];
    direction_pupil               = [];
    p_val_orig_vs_removed         = [];
    direction_orig_vs_removed     = [];
    
    if type_i > 2
        median_stationary_fr_ME    = [];
        median_motion_fr_ME        = [];
        p_val_ME                   = [];
        direction_ME               = [];
    end
    
    for probe_i = 1 : length(probe_ids)
        
        for clust_i = 1 : size(stationary_fr_store{probe_i}{type_i}, 2)
            
            stat = stationary_fr_store{probe_i}{type_i}(:, clust_i);
            mot = motion_fr_store{probe_i}{type_i}(:, clust_i);
            
            median_stationary_fr(end+1) = median(stat);
            median_motion_fr(end+1) = median(mot);
            [~, ~, p_val(end+1), direction(end+1)] = compare_groups_with_signrank(stat, mot);
            
            stat_pupil = stationary_fr_store_pupil{probe_i}{type_i}(:, clust_i);
            mot_pupil = motion_fr_store_pupil{probe_i}{type_i}(:, clust_i);
            
            median_stationary_fr_pupil(end+1) = median(stat_pupil);
            median_motion_fr_pupil(end+1) = median(mot_pupil);
            [~, ~, p_val_pupil(end+1), direction_pupil(end+1)] = compare_groups_with_signrank(stat_pupil, mot_pupil);
            
            
            if type_i > 2
                stat_ME = stationary_fr_store_ME{probe_i}{type_i}(:, clust_i);
                mot_ME = motion_fr_store_ME{probe_i}{type_i}(:, clust_i);
                
                median_stationary_fr_ME(end+1) = median(stat_ME);
                median_motion_fr_ME(end+1) = median(mot_ME);
                [~, ~, p_val_ME(end+1), direction_ME(end+1)] = compare_groups_with_signrank(stat_ME, mot_ME);
            end
            
            [~, ~, p_val_orig_vs_removed(end+1), direction_orig_vs_removed(end+1)] = compare_groups_with_signrank(mot, mot_pupil);
        end
        
        
    end
    
    if iscell(trial_types{type_i})
        type_str = trial_types{type_i}{1};
    else
        type_str = trial_types{type_i};
    end
    
    figure(1);
    h_ax = subplot(2, 2, type_i);
    hold on;
    
    fmt.xy_limits       = [0, 60];
    fmt.tick_space      = 20;
    fmt.line_order      = 'top';
    fmt.xlabel          = 'FR baseline';
    fmt.ylabel          = trial_types{type_i};
    fmt.include_inset   = false;
    fmt.colour_by       = 'significance';

    unity_plot_plot(h_ax, median_stationary_fr, median_motion_fr, direction, fmt);
    
    title(sprintf('%s, original', type_str), 'interpreter', 'none');
    
    
    figure(2);
    h_ax = subplot(2, 2, type_i);
    
    hold on;
    
    fmt.xy_limits       = [0, 60];
    fmt.tick_space      = 20;
    fmt.line_order      = 'top';
    fmt.xlabel          = 'FR baseline';
    fmt.ylabel          = trial_types{type_i};
    fmt.include_inset   = false;
    fmt.colour_by       = 'significance';

    unity_plot_plot(h_ax, median_stationary_fr_pupil, median_motion_fr_pupil, direction_pupil, fmt);
    
    title(sprintf('%s, saccades removed', type_str), 'interpreter', 'none');
    
    
    
    if type_i > 2
        
        figure(3);
        h_ax = subplot(2, 2, type_i);
        
        hold on;
        
        fmt.xy_limits       = [0, 60];
        fmt.tick_space      = 20;
        fmt.line_order      = 'top';
        fmt.xlabel          = 'FR baseline';
        fmt.ylabel          = trial_types{type_i};
        fmt.include_inset   = false;
        fmt.colour_by       = 'significance';
        
        unity_plot_plot(h_ax, median_stationary_fr_ME, median_motion_fr_ME, direction_ME, fmt);
        
        title(sprintf('%s, motion energy periods removed', type_str), 'interpreter', 'none');
    end
    
    
    
    figure(4);
    h_ax = subplot(2, 2, type_i);
    
    hold on;
    
    fmt.xy_limits       = [0, 60];
    fmt.tick_space      = 20;
    fmt.line_order      = 'top';
    fmt.xlabel          = sprintf('%s, original', type_str);
    fmt.ylabel          = sprintf('%s, removed', type_str);
    fmt.include_inset   = false;
    fmt.colour_by       = 'significance';

    unity_plot_plot(h_ax, median_motion_fr, median_motion_fr_pupil, direction_orig_vs_removed, fmt);
    
    title(sprintf('%s, saccades removed vs original', type_str), 'interpreter', 'none');
    
end


        
