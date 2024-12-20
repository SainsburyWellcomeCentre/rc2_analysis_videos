% body_motion_energy_and_saccade_analysis.m
%
% Description:
%   This script analyzes firing rate changes in neuron populations across
%   different experimental conditions by producing unity plots that reflect 
%   responses under three scenarios:
%       1. Original data, including all motion events.
%       2. Data with high-motion energy periods removed, particularly in replayed 
%          trials to isolate periods of low motion energy.
%       3. Data with periods just after pupil saccades removed, to examine 
%          post-saccadic effects on neural activity.
%
% Analysis Workflow:
%   BODY MOTION ENERGY
%       1. For each VT (visual task) and V (visual only) trial, the distribution 
%          of motion energy during running periods of the original trial (i.e., 
%          the trial from which the replay originates) is calculated.
%       2. A threshold motion value corresponding to a specific percentile 
%          is applied to body motion energy during both stationary and motion 
%          periods, retaining only time bins where motion energy falls below 
%          this threshold.
%       3. The firing rate (FR) is calculated from this subset of time bins 
%          and averaged across the trial for each cluster, then the median is 
%          calculated across trials.
%
%   SACCADE REMOVAL
%       1. For each session, a precomputed table of saccade frames is loaded.
%       2. For each trial, a time window is defined post-saccade using 
%          `duration_to_remove`. Time bins within this window are removed to 
%          create a saccade-exclusion mask.
%       3. The mask is applied to the firing rate data, and the mean across 
%          the trial and the median across trials are calculated for each cluster.
%
% Outputs:
%   - Unity plots comparing firing rates across stationary and motion periods 
%     in each scenario: original, high-motion energy excluded, and post-saccade 
%     periods removed.
%   - Summary statistics on data retention after removing post-saccade periods 
%     and high-motion energy segments.
%   - Modulation indices that quantify population responses in each condition.
%
% Usage:
%   Set the `experiment_groups` and `trial_types` variables as needed, and run 
%   the script to analyze firing rate distributions across conditions, comparing 
%   changes under different data exclusions.

close all;

experiment_groups           = 'visual_flow';
trial_types                 = {'RVT', 'RV', {'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};
duration_to_remove          = 0.25;     % duration after saccade to remove
motion_threshold_percentile = 5;

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
                
                % select the motion periods and get the lowest x-th prctile
                cam_motion_original         = original_motion_energy(original_motion_mask);
                motion_threshold            = prctile(cam_motion_original, motion_threshold_percentile);
                
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
%                 
%                 figure(trial_i);
%                 hold on;
%                 plot(original_motion_mask, 'r');
%                 plot(original_stationary_mask, 'b');
%                 plot(stationary_mask_ME * 1.2, 'y');
%                 plot(motion_mask_ME * 1.2, 'g');
%                 plot(replay_motion_energy, 'k');
%                 yline(motion_threshold);
%                 
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
    
    total_stationary_time = 0;
    total_motion_time = 0;
    
    total_stationary_saccade_time = 0;
    total_motion_saccade_time = 0;
    
    for type_i = 1 : length(trial_types)
        % sum all the the times of the trials together across types
        total_stationary_time = total_stationary_time + sum(stationary_time{probe_i}{type_i});
        total_motion_time = total_motion_time + sum(motion_time{probe_i}{type_i}(:, 1));
        
        total_stationary_saccade_time = total_stationary_saccade_time + sum(stationary_time_pupil{probe_i}{type_i}(:, 1));
        total_motion_saccade_time = total_motion_saccade_time + sum(motion_time_pupil{probe_i}{type_i}(:, 1));
    end
    
    sprintf('Probe: %s.\nData retained after saccades removal: %f (retained), %f (total), %f (percentage).', ...
        probe_ids{probe_i}, ...
        total_motion_saccade_time, ...
        total_motion_time, ...
        total_motion_saccade_time / total_motion_time)
    
    total_stationary_time = 0;
    total_motion_time = 0;
    total_stationary_ME_time = 0;
    total_motion_ME_time = 0;
    for type_i = 3 : length(trial_types)
        % only for VT and V
        total_stationary_time = total_stationary_time + sum(stationary_time{probe_i}{type_i});
        total_motion_time = total_motion_time + sum(motion_time{probe_i}{type_i}(:, 1));
        
        total_stationary_ME_time = total_stationary_ME_time + sum(stationary_time_ME{probe_i}{type_i}(:, 1));
        total_motion_ME_time = total_motion_ME_time + sum(motion_time_ME{probe_i}{type_i}(:, 1));
    end
    sprintf('\nData retained after body ME removal: %f (retained), %f (total), %f (percentage).', ...
        total_motion_ME_time, ...
        total_motion_time, ...
        total_motion_ME_time / total_motion_time)
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
    
    
end


%%

median_motion_fr_VT = [];
median_motion_fr_V = [];
direction = [];

median_motion_fr_VT_pupil = [];
median_motion_fr_V_pupil = [];
direction_pupil = [];

median_motion_fr_VT_ME = [];
median_motion_fr_V_ME = [];
direction_ME = [];

for probe_i = 1 : length(probe_ids)
        
        for clust_i = 1 : size(stationary_fr_store{probe_i}{type_i}, 2)
            

            % first plot, original
            mot_VT = motion_fr_store{probe_i}{3}(:, clust_i);
            mot_V = motion_fr_store{probe_i}{4}(:, clust_i);

            median_motion_fr_VT(end+1) = median(mot_VT);
            median_motion_fr_V(end+1) = median(mot_V);

            [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(mot_V, mot_VT);


            % first plot, saccade removal
            mot_VT_pupil = motion_fr_store_pupil{probe_i}{3}(:, clust_i);
            mot_V_pupil = motion_fr_store_pupil{probe_i}{4}(:, clust_i);

            median_motion_fr_VT_pupil(end+1) = median(mot_VT_pupil);
            median_motion_fr_V_pupil(end+1) = median(mot_V_pupil);

            [~, ~, ~, direction_pupil(end+1)] = compare_groups_with_signrank(mot_V_pupil, mot_VT_pupil);


            % first plot, saccade removal
            mot_VT_ME = motion_fr_store_ME{probe_i}{3}(:, clust_i);
            mot_V_ME = motion_fr_store_ME{probe_i}{4}(:, clust_i);

            median_motion_fr_VT_ME(end+1) = median(mot_VT_ME);
            median_motion_fr_V_ME(end+1) = median(mot_V_ME);

            [~, ~, ~, direction_ME(end+1)] = compare_groups_with_signrank(mot_V_ME, mot_VT_ME);

        end
end


figure(4);
h_ax = subplot(1, 3, 1);
hold on;

fmt.xy_limits       = [0, 60];
fmt.tick_space      = 20;
fmt.line_order      = 'top';
fmt.xlabel          = trial_types{4};
fmt.ylabel          = trial_types{3};
fmt.include_inset   = false;
fmt.colour_by       = 'significance';

unity_plot_plot(h_ax, median_motion_fr_V, median_motion_fr_VT, direction, fmt);

title('Original');


h_ax = subplot(1, 3, 2);
hold on;

fmt.xy_limits       = [0, 60];
fmt.tick_space      = 20;
fmt.line_order      = 'top';
fmt.xlabel          = trial_types{4};
fmt.ylabel          = trial_types{3};
fmt.include_inset   = false;
fmt.colour_by       = 'significance';

unity_plot_plot(h_ax, median_motion_fr_V_pupil, median_motion_fr_VT_pupil, direction_pupil, fmt);

title('No saccades');


h_ax = subplot(1, 3, 3);
hold on;

fmt.xy_limits       = [0, 60];
fmt.tick_space      = 20;
fmt.line_order      = 'top';
fmt.xlabel          = trial_types{4};
fmt.ylabel          = trial_types{3};
fmt.include_inset   = false;
fmt.colour_by       = 'significance';

unity_plot_plot(h_ax, median_motion_fr_V_ME, median_motion_fr_VT_ME, direction_ME, fmt);

title('No motion energy');


figure(5);

hold on;
modulation_index_no_saccades = [];
modulation_index_all_data = [];

for clust_i = 1 : 39
    modulation_index_no_saccades(end+1) = (median_motion_fr_VT_pupil(clust_i) - median_motion_fr_V_pupil(clust_i))...
        / (median_motion_fr_VT_pupil(clust_i) + median_motion_fr_V_pupil(clust_i));
    
    modulation_index_all_data(end+1) = (median_motion_fr_VT(clust_i) - median_motion_fr_V(clust_i))...
        / (median_motion_fr_VT(clust_i) + median_motion_fr_V(clust_i));
    
    if direction(clust_i) ~= 0
        if direction_pupil(clust_i) == 1
            scatter(2, modulation_index_no_saccades(clust_i), scatterball_size(3), 'red', 'o');
        elseif direction_pupil(clust_i) == -1
            scatter(2, modulation_index_no_saccades(clust_i), scatterball_size(3), 'blue', 'o');
        else 
            scatter(2, modulation_index_no_saccades(clust_i), scatterball_size(3), 'black', 'o');
        end

        if direction(clust_i) == 1
            scatter(1, modulation_index_all_data(clust_i), scatterball_size(3), 'red', 'o');  
        elseif direction(clust_i) == -1
            scatter(1, modulation_index_all_data(clust_i), scatterball_size(3), 'blue', 'o');
        end
        plot([1 2], [modulation_index_all_data(clust_i), modulation_index_no_saccades(clust_i)], 'black');
    end
end
xlim([0 3]);
ylim([-1.2 1.2]);
title('MI all data / no saccades');

only_responsive = direction ~= 0;
avg_mi_all_data = nanmean(modulation_index_all_data(only_responsive))
std_mi_all_data = nanstd(modulation_index_all_data(only_responsive))
sem_mi_all_data = nanstd(modulation_index_all_data(only_responsive)) / sqrt(39)

avg_mi_no_saccades  = nanmean(modulation_index_no_saccades(only_responsive))
std_mi_no_saccades  = nanstd(modulation_index_no_saccades(only_responsive))
sem_mi_no_saccades  = nanstd(modulation_index_no_saccades(only_responsive)) / sqrt(39)

[p] = signrank(modulation_index_all_data(only_responsive), modulation_index_no_saccades(only_responsive))


figure(6);

hold on;
modulation_index_no_body_ME = [];
modulation_index_all_data = [];

for clust_i = 1 : 39
    modulation_index_no_body_ME(end+1) = (median_motion_fr_VT_ME(clust_i) - median_motion_fr_V_ME(clust_i))...
        / (median_motion_fr_VT_ME(clust_i) + median_motion_fr_V_ME(clust_i));
    
    modulation_index_all_data(end+1) = (median_motion_fr_VT(clust_i) - median_motion_fr_V(clust_i))...
        / (median_motion_fr_VT(clust_i) + median_motion_fr_V(clust_i));
    
    if direction(clust_i) ~= 0
        if direction_ME(clust_i) == 1
            scatter(2, modulation_index_no_body_ME(clust_i), scatterball_size(3), 'red', 'o');
        elseif direction_ME(clust_i) == -1
            scatter(2, modulation_index_no_body_ME(clust_i), scatterball_size(3), 'blue', 'o');
        else 
            scatter(2, modulation_index_no_body_ME(clust_i), scatterball_size(3), 'black', 'o');
        end

        if direction(clust_i) == 1
            scatter(1, modulation_index_all_data(clust_i), scatterball_size(3), 'red', 'o');  
        elseif direction(clust_i) == -1
            scatter(1, modulation_index_all_data(clust_i), scatterball_size(3), 'blue', 'o');
        end
        plot([1 2], [modulation_index_all_data(clust_i), modulation_index_no_body_ME(clust_i)], 'black');
    end
end
xlim([0 3]);
ylim([-1.2 1.2]);
title('MI all data / no body ME');

only_responsive = direction ~= 0;
avg_mi_all_data = nanmean(modulation_index_all_data(only_responsive))
std_mi_all_data = nanstd(modulation_index_all_data(only_responsive))
sem_mi_all_data = nanstd(modulation_index_all_data(only_responsive))  / sqrt(39)

avg_mi_no_body_ME  = nanmean(modulation_index_no_body_ME(only_responsive))
std_mi_no_body_ME  = nanstd(modulation_index_no_body_ME(only_responsive))
sem_mi_no_body_ME  = nanstd(modulation_index_no_body_ME(only_responsive)) / sqrt(39)

[p] = signrank(modulation_index_all_data(only_responsive), modulation_index_no_body_ME(only_responsive))

















        
