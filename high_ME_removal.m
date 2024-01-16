%% Check firing rates during periods of high motion energy
%
%   compares firing rate during motion and stationary periods with and
%   without periods of high motion energy

%%
experiment_groups   = {'darkness', 'mismatch_darkness_oct21', 'mismatch_jul21'};
trial_types         = {{'T_RT', 'T_R', 'T'}, {'T_RT', 'T_R', 'T'}, {'T_RT', 'T_R', 'T'}};
fig_label           = 'T';
save_dir            = 'C:\Users\lee\Desktop\T';
save_on             = true;


%%
ctl                     = RC2Analysis();

prc_removed_motion      = {};
prc_removed_stationary  = {};
probe_id_store          = {};
experiment_group_store  = {};
noncropped_motion_fr    = {};
cropped_motion_fr       = {};
noncropped_stationary_fr = {};
cropped_stationary_fr   = {};
cluster_ids             = {};
spike_class             = {};
n_probe_recs            = 0;

for exp_i = 1 : length(experiment_groups)
    
    % probe recordings for this experiment group
    probe_ids = ctl.get_probe_ids(experiment_groups{exp_i});
    
    for ii = 1 : length(probe_ids)
        
        % skip this mouse as the camera was started late
        if strcmp(probe_ids{ii}, 'CAA-1114977_rec1_rec2_rec3')
            continue
        end
        
        % get probe recording data
        data                            = ctl.load_formatted_data(probe_ids{ii});
        clusters                        = data.VISp_clusters();
        
        % get trials of correct type
        trials                          = data.get_trials_with_trial_group_label(trial_types{exp_i});
        
        n_probe_recs                     = n_probe_recs + 1;
        
        % setup storage arrays
        prc_removed_motion{n_probe_recs}       = [];
        prc_removed_stationary{n_probe_recs}   = [];
        
        noncropped_motion_fr{n_probe_recs}     = [];
        cropped_motion_fr{n_probe_recs}        = [];
        
        noncropped_stationary_fr{n_probe_recs} = [];
        cropped_stationary_fr{n_probe_recs}    = [];
        
        cluster_ids{n_probe_recs}              = [];
        spike_class{n_probe_recs}              = [];
        
        probe_id_store{n_probe_recs}           = probe_ids{ii};
        experiment_group_store{n_probe_recs}   = experiment_groups{exp_i};
        
        for jj = 1 : length(trials)
            
            fprintf('%i/%i, %i/%i\n', jj, length(trials), ii, length(probe_ids));
            
            % get this trial and original trial
            replay_trial                = trials{jj}.to_aligned;
            original_trial              = replay_trial.original_trial;
            
            % mask of motion in the original trial
            original_motion_mask        = original_trial.motion_mask;
            original_stationary_mask    = original_trial.stationary_mask;
            
            % camera motion in original trial
            original_motion_energy      = original_trial.camera1;
            
            % select the motion periods and get the lowest 5th prctile
            cam_motion_original         = original_motion_energy(original_motion_mask);
            motion_threshold            = prctile(cam_motion_original, 5);
            
            % camera motion in replay trial, aligned
            replay_motion_energy        = replay_trial.camera1;
            
            % select the camera motion during treadmill motion periods in the replay trial
            cam_motion_replay           = replay_motion_energy(original_motion_mask);
            % and in stationary periods
            cam_stationary_replay       = replay_motion_energy(original_stationary_mask);
            
            % percentage of data removed due to the motion energy
            prc_removed_motion{n_probe_recs}(jj)     = get_prc_removed(cam_motion_replay > motion_threshold);
            prc_removed_stationary{n_probe_recs}(jj) = get_prc_removed(cam_stationary_replay > motion_threshold);
            
            % timebase of the replay trial
            replay_timebase                         = replay_trial.probe_t;
            replay_motion_mask                      = replay_trial.motion_mask;
            replay_stationary_mask                  = replay_trial.stationary_mask;
            
            % for each cluster
            for kk = 1 : length(clusters)
                
                % convolved firing rate during replay trial, aligned
                fr                                      = clusters(kk).fr.get_convolution(replay_timebase);
                
                % the original computed firing rate
                noncropped_motion_fr{n_probe_recs}(jj, kk)       = mean(fr(replay_motion_mask));
                noncropped_stationary_fr{n_probe_recs}(jj, kk)   = mean(fr(replay_stationary_mask));
                
                % the average firing rate with 'camera motion' removed
                cropped_motion_fr{n_probe_recs}(jj, kk)          = mean(fr(replay_motion_mask & replay_motion_energy < motion_threshold));
                cropped_stationary_fr{n_probe_recs}(jj, kk)      = mean(fr(replay_stationary_mask & replay_motion_energy < motion_threshold));
                
                % store info
                cluster_ids{n_probe_recs}(kk)                    = clusters(kk).id;
                spike_class{n_probe_recs}{kk}                    = clusters(kk).spiking_class;
            end
        end
    end
end



%%
print_statistics(noncropped_motion_fr, cropped_motion_fr, 'Motion');
print_statistics(noncropped_stationary_fr, cropped_stationary_fr, 'Stationary');


%%
% plot for each mouse and for each trial the fraction of data that will be removed

figure('position', [80, 180, 1500, 800]);

subplot(1, 2, 1);
hold on;
prc_removed_scatter(prc_removed_motion, probe_id_store, 'Motion');

subplot(1, 2, 2);
hold on;
prc_removed_scatter(prc_removed_stationary, probe_id_store, 'Stationary');

FigureTitle(gcf, fig_label);

format_and_save(gcf, save_on, save_dir, fig_label, 'prc_removed_per_mouse')



%% 
% pool fraction of data per trial that will be removed across mice
figure;
subplot(1, 2, 1);
prc_removed_histogram(prc_removed_motion, 'Motion');

subplot(1, 2, 2);
prc_removed_histogram(prc_removed_stationary, 'Stationary');

FigureTitle(gcf, fig_label);

format_and_save(gcf, save_on, save_dir, fig_label, 'prc_removed_histogram')




%%
% plot for each cluster the difference between firing during all "motion" periods and firing
% during "motion" periods with "motion frames" removed

motion_delta_fr             = cellfun(@(x, y)(x - y), noncropped_motion_fr, cropped_motion_fr, 'uniformoutput', false);
stationary_delta_fr         = cellfun(@(x, y)(x - y), noncropped_stationary_fr, cropped_stationary_fr, 'uniformoutput', false);

motion_delta_fr_median      = cellfun(@(x)(median(x)), motion_delta_fr, 'uniformoutput', false);
stationary_delta_fr_median  = cellfun(@(x)(median(x)), stationary_delta_fr, 'uniformoutput', false);

edges = -0.425:0.05:0.425;

figure;
subplot(1, 2, 1);
plot_delta_fr(motion_delta_fr_median, edges, 'Motion') 
xlabel(sprintf('\\Delta FR (Hz), (%s (all) - %s (cam motion removed))', fig_label, fig_label))
ylabel('# clusters');

subplot(1, 2, 2);
plot_delta_fr(stationary_delta_fr_median, edges, 'Stationary')

format_and_save(gcf, save_on, save_dir, fig_label, 'delta_fr')



%%
% plot for each cluster the difference between firing during all "motion"
% periods and firing during "motion" periods with "motion frames" removed,
% but restricted to trials in which at least some data has been removed

idx_non_zero_motion     = cellfun(@(x)(x > 0), prc_removed_motion, 'uniformoutput', false);
idx_non_zero_stationary = cellfun(@(x)(x > 0), prc_removed_stationary, 'uniformoutput', false);

motion_delta_fr         = cellfun(@(x, y, z)(x(z, :) - y(z, :)), ...
    noncropped_motion_fr, cropped_motion_fr, idx_non_zero_motion, 'uniformoutput', false);
stationary_delta_fr     = cellfun(@(x, y, z)(x(z, :) - y(z, :)), ...
    noncropped_stationary_fr, cropped_stationary_fr, idx_non_zero_stationary, 'uniformoutput', false);

motion_delta_fr_median      = cellfun(@(x)(median(x)), motion_delta_fr, 'uniformoutput', false);
stationary_delta_fr_median  = cellfun(@(x)(median(x)), stationary_delta_fr, 'uniformoutput', false);


figure;
subplot(1, 2, 1);
plot_delta_fr(motion_delta_fr_median, edges, 'Motion') 
xlabel(sprintf('\\Delta FR (Hz), (%s (all) - %s (cam motion removed))', fig_label, fig_label))
ylabel('# clusters');

subplot(1, 2, 2);
plot_delta_fr(stationary_delta_fr_median, edges, 'Stationary')

format_and_save(gcf, save_on, save_dir, fig_label, 'delta_fr_trials_gt_zero_prc')



%% AUX FUNCTIONS

function prc = get_prc_removed(mask)

n_samples_above = sum(mask);
n_samples_total = length(mask);
prc = 100 * n_samples_above / n_samples_total;

end



function prc_removed_scatter(prc_removed, probe_id_store, title_str)

for ii = 1 : length(prc_removed)
    scatter(ii, prc_removed{ii}, [], 'k', 'fill');    
end

xlim([0, length(prc_removed)+1])
ylim([0, 100])

set(gca, 'xtick', 1:length(prc_removed), ...
         'xticklabel', probe_id_store, ...
         'ticklabelinterpreter', 'none', ...
         'plotboxaspectratio', [3, 1, 1]);

xtickangle(30);

title(title_str);

end



function prc_removed_histogram(prc_removed, title_str)

histogram([prc_removed{:}], 'binwidth', 5);
M = max(60, max([prc_removed{:}]));
set(gca, 'xlim', [0, M], 'ylim', [0, 100]);
box off;
xlabel('% motion period removed');
ylabel('# trials');
title(title_str);

end



function plot_delta_fr(delta_fr, edges, title_str)

histogram([delta_fr{:}], 'binedges', edges);

p_val = signtest([delta_fr{:}]);

text(0.5, 50, sprintf('p = %.2f', p_val), 'HorizontalAlignment', 'right', 'verticalalignment', 'top');
box off;
title(title_str);

end



function format_and_save(h_fig, save_on, save_dir, fig_label, suffix)

set(h_fig,  'paperunits', 'inches', ...
    'papersize', [15, 10], ...
    'paperposition', [0, 0, 15, 10], ...
    'renderer', 'painters');

if save_on
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    fname = fullfile(save_dir, sprintf('%s_%s.pdf', fig_label, suffix));
    print(fname, '-dpdf');
end

end



function print_statistics(noncropped_fr, cropped_fr, str)

delta_fr            = [];
noncropped_median  = [];
cropped_median      = [];

for probe_i = 1 : length(noncropped_fr)
    
    trial_delta_fr = noncropped_fr{probe_i} - cropped_fr{probe_i};
    delta_fr = [delta_fr, median(trial_delta_fr, 1)];
    
    noncropped_median = [noncropped_median, median(noncropped_fr{probe_i}, 1)];
    cropped_median = [cropped_median, median(cropped_fr{probe_i}, 1)];
end

p_val_signrank = signrank(noncropped_median(:), cropped_median(:));
p_val_sign = signtest(delta_fr);

fprintf('%s\n', str);
fprintf('p-value signrank: %.2f\n', p_val_signrank);
fprintf('p-value sign-test: %.2f\n', p_val_sign);
fprintf('\n\n\n');
end
