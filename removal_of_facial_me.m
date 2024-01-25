experiment_groups           = 'visual_flow';
trial_types                 = {{'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};
duration_to_remove          = 0.25;   
motion_threshold_percentile = 5;

ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);


median_motion_V  = [];
median_motion_VT  = [];
direction = [];

for probe_i = 1 : length(probe_ids)
    data   = ctl.load_formatted_data(probe_ids{probe_i});
    clusters  = data.VISp_clusters();
    
    % =====================================================================
    % analyse VT
    type_i = 1; 
    
    % Get the distribution of facial ME in stationary and motion periods
    trials = data.get_trials_with_trial_group_label(trial_types{type_i});
    facial_me_motion_all = zeros(length(trials), 200000);
    facial_me_motion_all_stationary = zeros(length(trials), 200000);
    for trial_i = 1 : length(trials)
        trial  = trials{trial_i}.to_aligned;
        original_trial              = trial.original_trial;

        original_motion_mask        = original_trial.motion_mask;
        original_stationary_mask    = original_trial.stationary_mask;

        face_motion_energy      = trial.camera0;
        face_motion_energy_masked  = face_motion_energy(original_motion_mask);
        face_motion_energy_masked_stationary  = face_motion_energy(original_stationary_mask);

        facial_me_motion_all(trial_i, 1:length(face_motion_energy_masked)) = face_motion_energy_masked;
        facial_me_motion_all_stationary(trial_i, 1:length(face_motion_energy_masked_stationary)) = face_motion_energy_masked_stationary;
% 
%             figure(trial_i);
%             hold on;
%             histogram(cam_motion_original_stationary);
%             histogram(cam_motion_original_motion);
    end
    
    % Set the threshold
    facial_me_motion_all(facial_me_motion_all==0) = NaN;
    facial_me_motion_all_stationary(facial_me_motion_all_stationary==0) = NaN;
    no_facial_movements_threshold  = prctile(facial_me_motion_all(:), motion_threshold_percentile);


%         figure(type_i);
%         hold on;
%         histogram(facial_me_motion_all(:));
%         hold on;
%         histogram(facial_me_motion_all_stationary(:))
%         xline(no_facial_movements_threshold)
    
    % Calculate windows in which facial ME is low and the animal is running
    % Get the mean firing rate
    % Save the windows in a variables to be reused to analyse V
    windows_fme = zeros(length(trials), 350000);
    mean_spikes_VT = zeros(length(trials), length(clusters));
    for trial_i = 1 : length(trials)
        trial  = trials{trial_i}.to_aligned;
        original_trial              = trial.original_trial;

        original_motion_mask        = original_trial.motion_mask;

        face_motion_energy      = trial.camera0;
        f_me_mask = face_motion_energy < no_facial_movements_threshold;
        f_me_doubled_masking = f_me_mask & original_motion_mask(1:length(f_me_mask));
        windows_fme(trial_i, 1:length(f_me_mask)) = f_me_doubled_masking;

%         figure(trial_i + 1);
%         hold on;
%         plot(face_motion_energy);
%         plot(f_me_mask * 1.5);
%         plot(original_motion_mask * 2);
%         plot(f_me_doubled_masking);

        for clust_i = 1 : length(clusters)     
            fr = clusters(clust_i).fr.get_convolution(trial.probe_t);
            mean_spikes_VT(trial_i, clust_i) = mean(fr(f_me_doubled_masking));
        end
    end
    
    
    % Analyse V
    type_i = 2; 
    trials = data.get_trials_with_trial_group_label(trial_types{type_i});
    mean_spikes_V = zeros(length(trials), length(clusters));
    for trial_i = 1 : length(trials)
        trial  = trials{trial_i}.to_aligned;
        
        for clust_i = 1 : length(clusters)     
            fr = clusters(clust_i).fr.get_convolution(trial.probe_t);
            mean_spikes_V(trial_i, clust_i) = mean(fr(logical(windows_fme(trial_i, 1:length(fr)))));
        end
    end
    for clust_i = 1 : length(clusters) 
        mot_V = mean_spikes_V(:, clust_i);
        mot_VT = mean_spikes_VT(:, clust_i);
        median_motion_V(end+1) = median(mot_V);
        median_motion_VT(end+1) = median(mot_VT);
        [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(mot_V, mot_VT);
    end
end

figure(1);
h_ax = subplot(1, 1, 1);
hold on;
fmt.xy_limits       = [0, 60];
fmt.tick_space      = 20;
fmt.line_order      = 'top';
fmt.xlabel          = trial_types{2};
fmt.ylabel          = trial_types{2};
fmt.include_inset   = false;
fmt.colour_by       = 'significance';

unity_plot_plot(h_ax, median_motion_V, median_motion_VT, direction, fmt);


















