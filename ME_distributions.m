% trial_types                 = {'RVT_gain_up', 'RV_gain_up'}; %{'RVT', 'RV'};
trial_types                 = {{'RVT', 'RVT_gain_up'}, {'RV', 'RV_gain_up'}};

ctl                         = RC2Analysis();
% probe_ids                   = ctl.get_probe_ids('mismatch_nov20', 'mismatch_jul21');
probe_ids                   = ctl.get_probe_ids('visual_flow', 'mismatch_nov20', 'mismatch_jul21');

figure(1)
for probe_i = 1 : length(probe_ids)
    data                    = ctl.load_formatted_data(probe_ids{probe_i});

    distribution_RVT = [];
    trials = data.get_trials_with_trial_group_label(trial_types{1});
    for trial_i = 1 : length(trials)
        trial = trials{trial_i}.to_aligned;
        distribution_RVT = cat(1, distribution_RVT, trial.camera1(trial.motion_mask));
    end

    distribution_RV = [];
    trials = data.get_trials_with_trial_group_label(trial_types{2});
    for trial_i = 1 : length(trials)
        trial = trials{trial_i}.to_aligned;
        distribution_RV = cat(1, distribution_RV, trial.camera1(trial.motion_mask));
    end

    subplot(length(probe_ids), 1, probe_i)
    
    histogram(distribution_RVT)
    hold on
    histogram(distribution_RV)
    legend('RVT', 'RV')
    title(sprintf(probe_ids{probe_i}))
    
    sprintf('Probe: %s, Median RVT: %.2f', probe_ids{probe_i}, median(distribution_RVT))
    sprintf('Median RV: %.2f', median(distribution_RV))   
end