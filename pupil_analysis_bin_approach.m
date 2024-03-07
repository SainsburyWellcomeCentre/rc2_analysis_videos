close all;

experiment_groups           = 'visual_flow';
trial_types                 = {{'VT_RVT', 'VT_RV'}, {'V_RVT', 'V_RV'}};

% Thresholds for unity plot, make masks on VF + T
% low_percentile = 10;
% hig_percentile = 90;
custom_thresholds = [41.9419 61.4615 73.6737 57.5576];

% Thresholds for the "bins plot"
bin_num = 4;
bin_size = 100 / bin_num;


ctl                         = RC2Analysis();
probe_ids                   = ctl.get_probe_ids(experiment_groups);


median_VT = [];
median_V = [];
direction = []; 

median_VT_no_mask = [];
median_V_no_mask = [];
direction_no_mask = []; 

VT_fr_per_bin = zeros(39, 4);
V_fr_per_bin = zeros(39, 4);
clust_idx = 1;
for probe_i = 1 : length(probe_ids)
    data   = ctl.load_formatted_data(probe_ids{probe_i});
    clusters  = data.VISp_clusters();
    
    % =====================================================================
    pupil_diameter_motion_all = zeros(length(trial_types), 100, 300000);
    for type_i = 1 : length(trial_types)
        % Get the distribution of pupil diameter
        trials = data.get_trials_with_trial_group_label(trial_types{type_i});
        
        for trial_i = 1 : length(trials)
            trial                       = trials{trial_i}.to_aligned;
            original_trial              = trial.original_trial;
            
            original_motion_mask        = original_trial.motion_mask;

            pupil_diameter              = trial.pupil_diameter;
            pupil_diameter_masked       = pupil_diameter(original_motion_mask);

            pupil_diameter_motion_all(type_i, trial_i, 1:length(pupil_diameter_masked)) = pupil_diameter_masked;
        end
    end
    
    % Set the threshold
    pupil_diameter_motion_all(pupil_diameter_motion_all==0) = NaN;
%     lower_threshold  = prctile(pupil_diameter_motion_all(2, :), low_percentile);
%     high_threshold  = prctile(pupil_diameter_motion_all(2, :), hig_percentile);
    mask_threshold = custom_thresholds(probe_i);
    
    % four parts thresholds for plotting
    thresholds = zeros(length(trial_types), bin_num);
    for type_i = 1 : length(trial_types)
        for bin = 1 : bin_num
            thresholds(type_i, bin) = prctile(pupil_diameter_motion_all(type_i, :), bin_size * bin);
        end
    end
    
    edges = linspace(0, 100, 1000);
    figure(probe_i);
    hold on;
    title(probe_ids(probe_i));
    subplot(1, 2, 1);
    histogram(pupil_diameter_motion_all(1, :), edges);
    hold on;
    histogram(pupil_diameter_motion_all(2, :), edges);
    xlabel('Pupil diameter (pixel)');
    ylabel('Counts');
    xline(mask_threshold)
     
    mean_spikes = zeros(length(trial_types), length(trials), length(clusters));
	mean_spikes_no_mask = zeros(length(trial_types), length(trials), length(clusters));

    pd_doubled_masking = zeros(length(trials), 300000);
    
    fr_per_diameter = NaN(length(trial_types), length(clusters), bin_num, 300000);
    fr_per_diameter_mean = zeros(length(trial_types), length(clusters), bin_num);

    for type_i = 1 : length(trial_types)
        % Get the distribution of pupil diameter
        trials = data.get_trials_with_trial_group_label(trial_types{type_i});
        
        for trial_i = 1 : length(trials)
            trial  = trials{trial_i}.to_aligned;
            original_trial              = trial.original_trial;
            original_motion_mask        = original_trial.motion_mask;
            pupil_diameter              = trial.pupil_diameter;
%             pd_mask = (pupil_diameter < high_threshold) & (pupil_diameter > lower_threshold);
            pd_mask = pupil_diameter < mask_threshold;
            
            if type_i == 1
                % calcumate motion + range pupil in VF + T and use it in VF
                pd_doubled_masking(trial_i, 1:length(pd_mask)) = pd_mask & original_motion_mask(1:length(pd_mask));
            end
            
            for clust_i = 1 : length(clusters)     
                fr = clusters(clust_i).fr.get_convolution(trial.probe_t);
                mask = logical(pd_doubled_masking(trial_i, 1:length(fr)));                
                mean_spikes(type_i, trial_i, clust_i) = nanmean(fr(mask));
                mean_spikes_no_mask(type_i, trial_i, clust_i) = nanmean(fr(original_motion_mask));
                
%                 pupil_diameter_masked = pupil_diameter(original_motion_mask);
%                 fr_masked = fr(original_motion_mask);
%                 for point = 1 : length(fr_masked)
%                    if ~isnan(pupil_diameter_masked(point))
%                        bins = find((thresholds(type_i, :) >= pupil_diameter_masked(point)));
%                        fr_per_diameter(type_i, clust_i, bins(1), point) = fr_masked(point);
%                    end
%                 end
            end
        end
    end
    
    for clust_i = 1 : length(clusters)
        mean_VT = mean_spikes(1, :, clust_i);
        mean_V = mean_spikes(2, :, clust_i);
        median_VT(end+1) = nanmedian(mean_VT);
        median_V(end+1) = nanmedian(mean_V);
        [~, ~, ~, direction(end+1)] = compare_groups_with_signrank(mean_V, mean_VT);
        
        mean_VT_no_mask = mean_spikes_no_mask(1, :, clust_i);
        mean_V_no_mask = mean_spikes_no_mask(2, :, clust_i);
        median_VT_no_mask(end+1) = nanmedian(mean_VT_no_mask);
        median_V_no_mask(end+1) = nanmedian(mean_V_no_mask);
        [~, ~, ~, direction_no_mask(end+1)] = compare_groups_with_signrank(mean_V_no_mask, mean_VT_no_mask);
        
%         
%         for type_i = 1 : length(trial_types)
%             for bin = 1 : bin_num
%                 fr_per_diameter_mean(type_i, clust_i, bin) = nanmean(fr_per_diameter(type_i, clust_i, bin, :));
%             end
%         end
    end
%     subplot(1, 2, 2);
%     hold on;
%     for bin = 1 : bin_num
%         scatter(thresholds(1, bin), fr_per_diameter_mean(1, :, bin), 'blue');
%         scatter(thresholds(2, bin), fr_per_diameter_mean(2, :, bin), 'red');
%     end
%     
%     for clust_i = 1 : length(clusters)
%         for bin = 1 : bin_num
%             VT_fr_per_bin(clust_idx, bin) = fr_per_diameter_mean(1, clust_i, bin);
%             V_fr_per_bin(clust_idx, bin) = fr_per_diameter_mean(2, clust_i, bin);
%         end
%         clust_idx = clust_idx + 1;
%     end
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

unity_plot_plot(h_ax, median_V, median_VT, direction, fmt);


% 
% figure(6);
% hold on;
% for bin = 1 : bin_num
%     errorbar(bin, mean(VT_fr_per_bin(:, bin)), std(VT_fr_per_bin(:, bin)),'-bo','MarkerSize',10,...
%     'MarkerEdgeColor','b','MarkerFaceColor','b');
%     errorbar(bin, mean(V_fr_per_bin(:, bin)), std(V_fr_per_bin(:, bin)), '-ro','MarkerSize',10,...
%     'MarkerEdgeColor','r','MarkerFaceColor','r');
% end
% xlim([0 5]);


figure(7);
hold on;
modulation_index = [];
modulation_index_no_mask = [];
for clust_i = 1 : 39
    modulation_index(end+1) = (median_VT(clust_i) - median_V(clust_i)) / (median_VT(clust_i) + median_V(clust_i));
    modulation_index_no_mask(end+1) = (median_VT_no_mask(clust_i) - median_V_no_mask(clust_i)) / (median_VT_no_mask(clust_i) + median_V_no_mask(clust_i));
    
    if direction_no_mask(clust_i) ~= 0
        if direction(clust_i) == 1
            scatter(2, modulation_index(clust_i), scatterball_size(3), 'red', 'o');
        elseif direction(clust_i) == -1
            scatter(2, modulation_index(clust_i), scatterball_size(3), 'blue', 'o');
        else 
            scatter(2, modulation_index(clust_i), scatterball_size(3), 'black', 'o');
        end

        if direction_no_mask(clust_i) == 1
            scatter(1, modulation_index_no_mask(clust_i), scatterball_size(3), 'red', 'o');  
        elseif direction_no_mask(clust_i) == -1
            scatter(1, modulation_index_no_mask(clust_i), scatterball_size(3), 'blue', 'o');
        end
        plot([1 2], [modulation_index_no_mask(clust_i), modulation_index(clust_i)], 'black');
    end
end
xlim([0 3]);
ylim([-1.2 1.2]);


only_responsive_no_mask = direction_no_mask ~= 0;
avg_mi_no_mask = nanmean(modulation_index_no_mask(only_responsive_no_mask))
std_mi_no_mask = nanstd(modulation_index_no_mask(only_responsive_no_mask))
avg_mi  = nanmean(modulation_index(only_responsive_no_mask))
std_mi  = nanstd(modulation_index(only_responsive_no_mask))
[p] = signrank(modulation_index_no_mask(only_responsive_no_mask), modulation_index(only_responsive_no_mask))


[p_VT,tbl_VT,stats_VT] = anova1(VT_fr_per_bin);
[p_V,tbl_V,stats_V] = anova1(V_fr_per_bin);
