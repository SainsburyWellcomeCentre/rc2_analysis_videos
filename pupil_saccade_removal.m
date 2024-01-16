%% Check firing rates just after saccade removal
%
%   compares firing rate during motion and stationary periods with and
%   without the period just after saccade


%%
experiment_groups   = 'visual_flow';
trial_types         = {'RVT', 'RV', {'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};
duration_to_remove  = 0.25;  % amount of time in seconds to removal after a saccade


%%
ctl                     = RC2Analysis();

% probe recordings for this experiment group
probe_ids               = ctl.get_probe_ids(experiment_groups);

probe_id_store          = cell(1, length(trial_types));

prc_removed_motion      = cell(1, length(trial_types));
prc_removed_stationary  = cell(1, length(trial_types));
        
noncropped_motion_fr    = cell(1, length(trial_types));
cropped_motion_fr       = cell(1, length(trial_types));

noncropped_stationary_fr = cell(1, length(trial_types));
cropped_stationary_fr    = cell(1, length(trial_types));

cluster_ids             = cell(1, length(trial_types));
spike_class             = cell(1, length(trial_types));


for probe_i = 1 : length(probe_ids)
    
    % skip this mouse as the camera was started late
    if strcmp(probe_ids{probe_i}, 'CAA-1114977_rec1_rec2_rec3')
        continue
    end
    
    % get probe recording data
    data                            = ctl.load_formatted_data(probe_ids{probe_i});
    clusters                        = data.VISp_clusters();
    
    % get table of saccades
    sessions                        = data.motion_sessions;
    saccade_tbl                     = ctl.load.camera0_saccades(sessions{1}.session_id);
    
    % time of each saccade
    saccade_times                   = sessions{1}.camera_t(saccade_tbl.saccade_frame);
    
    for type_i = 1 : length(trial_types)
        
        % get trials of correct type
        trials                              = data.get_trials_with_trial_group_label(trial_types{type_i});
        
        % setup storage arrays
        prc_removed_motion{type_i}{probe_i}       = nan(length(trials), 1);
        prc_removed_stationary{type_i}{probe_i}   = nan(length(trials), 1);
        
        noncropped_motion_fr{type_i}{probe_i}     = nan(length(trials), length(clusters));
        cropped_motion_fr{type_i}{probe_i}        = nan(length(trials), length(clusters));
        
        noncropped_stationary_fr{type_i}{probe_i} = nan(length(trials), length(clusters));
        cropped_stationary_fr{type_i}{probe_i}    = nan(length(trials), length(clusters));
        
        cluster_ids{type_i}{probe_i}              = nan(length(clusters), 1);
        spike_class{type_i}{probe_i}              = cell(length(clusters), 1);
        
        probe_id_store{type_i}{probe_i}           = probe_ids{probe_i};
        
        figure
        
        for trial_i = 1 : length(trials)
            
            fprintf('%i/%i, %i/%i\n', trial_i, length(trials), probe_i, length(probe_ids));
            
            % get this trial and original trial
            trial                       = trials{trial_i}.to_aligned;
            
            % mask of motion in the trial
            motion_mask        = trial.motion_mask;
            stationary_mask    = trial.stationary_mask;
            
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
            
            % plot the saccades on the trials
            subplot(4, 5, trial_i)
            plot(trial.velocity)
            hold on;
            plot(saccade_mask * 5)
            box off
            if iscell(trial_types{type_i})
                str = trial_types{type_i}{1};
            else
                str = trial_types{type_i};
            end
            title(sprintf('%s, trial ID %i, %s', probe_ids{probe_i}, trial.trial_id, str), 'interpreter', 'none');
            
            % percentage of data removed due to the motion energy
            prc_removed_motion{type_i}{probe_i}(trial_i)     = sum(motion_mask & saccade_mask) / sum(motion_mask);
            prc_removed_stationary{type_i}{probe_i}(trial_i) = sum(stationary_mask & saccade_mask) / sum(stationary_mask);
             
            % for each cluster
            for clust_i = 1 : length(clusters)
                
                % convolved firing rate during replay trial, aligned
                fr                                              = clusters(clust_i).fr.get_convolution(timebase);
                
                % the original computed firing rate
                noncropped_motion_fr{type_i}{probe_i}(trial_i, clust_i)       = mean(fr(motion_mask));
                noncropped_stationary_fr{type_i}{probe_i}(trial_i, clust_i)   = mean(fr(stationary_mask));
                
                % the average firing rate with 'camera motion' removed
                cropped_motion_fr{type_i}{probe_i}(trial_i, clust_i)          = mean(fr(motion_mask & ~saccade_mask));
                cropped_stationary_fr{type_i}{probe_i}(trial_i, clust_i)      = mean(fr(stationary_mask & ~saccade_mask));
                
                % store info
                cluster_ids{type_i}{probe_i}(clust_i)                    = clusters(clust_i).id;
                spike_class{type_i}{probe_i}{clust_i}                    = clusters(clust_i).spiking_class;
            end
        end
    end
end 



%% print statistics
print_statistics(noncropped_motion_fr, cropped_motion_fr, trial_types, 'Motion')
print_statistics(noncropped_stationary_fr, cropped_stationary_fr, trial_types, 'Stationary')



function print_statistics(noncropped_fr, cropped_fr, trial_types, str)

delta_fr            = cell(1, length(trial_types));
noncropped_median  = cell(1, length(trial_types));
cropped_median      = cell(1, length(trial_types));

for type_i = 1 : length(trial_types)
    for probe_i = 1 : length(noncropped_fr)
    
        trial_delta_fr = noncropped_fr{type_i}{probe_i} - cropped_fr{type_i}{probe_i};
        delta_fr{type_i} = [delta_fr{type_i}, median(trial_delta_fr, 1)];
        
        noncropped_median{type_i} = [noncropped_median{type_i}, median(noncropped_fr{type_i}{probe_i}, 1)];
        cropped_median{type_i} = [cropped_median{type_i}, median(cropped_fr{type_i}{probe_i}, 1)];
    end
end

for type_i = 1 : length(trial_types)
    
    p_val_signrank = signrank(noncropped_median{type_i}(:), cropped_median{type_i});
    p_val_sign = signtest(delta_fr{type_i});
    
    % print info
    if iscell(trial_types{type_i})
        trial_str = strjoin(trial_types{type_i}, ',');
    else
        trial_str = trial_types{type_i};
    end
    
    fprintf('Period: %s\n', str);
    fprintf('Trial type: %s\n', trial_str);
    fprintf('  p-value Wilcoxon signrank: %.2f\n', p_val_signrank);
    fprintf('  p-value sign-test: %.2f\n', p_val_sign);
end

end