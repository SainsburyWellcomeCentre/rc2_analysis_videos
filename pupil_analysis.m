% How is this analysis done?
% 1. Get the distribution of pupil diameter for all VT trials in stationary
% and motion period (separated)
% 2. calculate the pupil diameter threshold which retains a percentile of
% stationary period data
% 3. on a trial by trial basis, calculate a pupil mask such that the pupil
% diameter is lower than the above threshold. 
% 4. we combine the pupil diameter mask with the motion mask (double mask) and we keep
% data under these two.
% 5. we apply the double mask to the V traces
% 6. we compare the firing rates with the same double mask (mean across
% trial per cluser, median over trials per cluster) 





experiment_groups           = 'visual_flow';
trial_types                 = {{'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};
percentile = 70;

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);


median_pd_V  = [];
median_pd_VT  = [];
direction = [];

for probe_i = 1 : length(probe_ids)
    data   = ctl.load_formatted_data(probe_ids{probe_i});
    clusters  = data.VISp_clusters();
    
    % =====================================================================
    % analyse VT
    type_i = 1; 
    
    % Get the distribution of pupil diameter in stationary and motion periods
    trials = data.get_trials_with_trial_group_label(trial_types{type_i});
    pupil_diameter_motion_all = zeros(length(trials), 200000);
    pupil_diameter_motion_all_stationary = zeros(length(trials), 200000);
    for trial_i = 1 : length(trials)
        trial  = trials{trial_i}.to_aligned;
        original_trial              = trial.original_trial;

        original_motion_mask        = original_trial.motion_mask;
        original_stationary_mask    = original_trial.stationary_mask;

        pupil_diameter      = trial.pupil_diameter;
        pupil_diameter_masked  = pupil_diameter(original_motion_mask);
        pupil_diameter_masked_stationary  = pupil_diameter(original_stationary_mask);

        pupil_diameter_motion_all(trial_i, 1:length(pupil_diameter_masked)) = pupil_diameter_masked;
        pupil_diameter_motion_all_stationary(trial_i, 1:length(pupil_diameter_masked_stationary)) = pupil_diameter_masked_stationary;

        %     distribution per trial plot
%             figure(trial_i);
%             hold on;
%             histogram(cam_motion_original_stationary);
%             histogram(cam_motion_original_motion);
    end
    
    % Set the threshold
    pupil_diameter_motion_all(pupil_diameter_motion_all==0) = NaN;
    pupil_diameter_motion_all_stationary(pupil_diameter_motion_all_stationary==0) = NaN;
    small_diameter_threshold  = prctile(pupil_diameter_motion_all_stationary(:), percentile);


    figure(probe_i);
    hold on;
    h = histogram(pupil_diameter_motion_all(:), 50);
    hold on;
    g = histogram(pupil_diameter_motion_all_stationary(:), 50);
    xline(small_diameter_threshold)
    
    % Save the windows in a variables to be reused to analyse V
    windows_pd = zeros(length(trials), 350000);
    mean_spikes_VT = zeros(length(trials), length(clusters));
    
    total_pupil_diameter_mask_len = 0;
    total_double_mask_len = 0;
    for trial_i = 1 : length(trials)
        trial  = trials{trial_i}.to_aligned;
        original_trial              = trial.original_trial;

        original_motion_mask        = original_trial.motion_mask;

        pupil_diameter      = trial.pupil_diameter;
        pd_mask = pupil_diameter < small_diameter_threshold;
        
        pd_doubled_masking = pd_mask & original_motion_mask(1:length(pd_mask));
        motion_mask = original_motion_mask(1:length(pd_mask));
        windows_pd(trial_i, 1:length(pd_mask)) = pd_doubled_masking;
        
%         sprintf('Total data length: %f; retained: %f; discarted: %f', ...
%                                     length(pupil_diameter(motion_mask)), ...
%                                     length(pupil_diameter(pd_doubled_masking)) / length(pupil_diameter(motion_mask)), ...
%                                     1 - length(pupil_diameter(pd_doubled_masking)) / length(pupil_diameter(motion_mask))) 
                                
        total_pupil_diameter_mask_len = total_pupil_diameter_mask_len + length(pupil_diameter(motion_mask));
        total_double_mask_len = total_double_mask_len + length(pupil_diameter(pd_doubled_masking));
        
        % To see masks for each trial, uncomment this section                    
%         figure(trial_i);
%         hold on;
%         plot(pupil_diameter/70, 'g');
%         plot(pd_mask * 1, 'k');
%         plot(original_motion_mask * 2, 'r');
%         plot(pd_doubled_masking, 'y');
%         plot(original_trial.stationary_mask * 2, 'b')


        for clust_i = 1 : length(clusters)     
            fr = clusters(clust_i).fr.get_convolution(trial.probe_t);
            mean_spikes_VT(trial_i, clust_i) = nanmean(fr(pd_doubled_masking));
        end
    end
    
    
    sprintf('Probe %s; percentage of data retained with the pupil masking: %f; threshold: %f', ...
        probe_ids{probe_i}, total_double_mask_len / total_pupil_diameter_mask_len, small_diameter_threshold)
    
    % =====================================================================
    % Analyse V
    type_i = 2; 
    trials = data.get_trials_with_trial_group_label(trial_types{type_i});
    mean_spikes_V = zeros(length(trials), length(clusters));
    for trial_i = 1 : length(trials)
        trial  = trials{trial_i}.to_aligned;
        
        for clust_i = 1 : length(clusters)     
            fr = clusters(clust_i).fr.get_convolution(trial.probe_t);
            mean_spikes_V(trial_i, clust_i) = nanmean(fr(logical(windows_pd(trial_i, 1:length(fr)))));
        end
    end
    for clust_i = 1 : length(clusters) 
        pd_V = mean_spikes_V(:, clust_i);
        pd_VT = mean_spikes_VT(:, clust_i);
        median_pd_V(end+1) = nanmedian(pd_V);
        median_pd_VT(end+1) = nanmedian(pd_VT);
        [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(pd_V, pd_VT);
    end
end

figure(5);
h_ax = subplot(1, 1, 1);
hold on;
fmt.xy_limits       = [0, 60];
fmt.tick_space      = 20;
fmt.line_order      = 'top';
fmt.xlabel          = trial_types{2};
fmt.ylabel          = trial_types{1};
fmt.include_inset   = false;
fmt.colour_by       = 'significance';

unity_plot_plot(h_ax, median_pd_V, median_pd_VT, direction, fmt);


















