function [mask_threshold, smooth_counts1, smooth_counts2] = find_mask_threshold(motion, stationary, edges, smooth_bins)
    % find_mask_threshold
    %
    % Description:
    %   This function calculates a threshold value (`mask_threshold`) by analyzing
    %   two distributions (motion and stationary data). It smooths the histogram 
    %   counts of each distribution and finds the intersection point after the 
    %   peak of the stationary distribution. This threshold helps in differentiating
    %   the two distributions for masking purposes.
    %
    % Inputs:
    %   - motion: Array of data points representing the motion distribution.
    %   - stationary: Array of data points representing the stationary distribution.
    %   - edges: Bin edges used to calculate histograms for both distributions.
    %   - smooth_bins: Number of bins used for smoothing the histograms.
    %
    % Outputs:
    %   - mask_threshold: The threshold value for masking based on the intersection
    %     point between the two smoothed histograms.
    %   - smooth_counts1: Smoothed histogram counts of the motion distribution.
    %   - smooth_counts2: Smoothed histogram counts of the stationary distribution.
    %
    % Usage:
    %   [mask_threshold, smooth_counts1, smooth_counts2] = find_mask_threshold(motion, stationary, edges, smooth_bins);
    
        % Compute histograms for motion and stationary data
        counts1 = histcounts(motion, edges);       % Histogram for motion data
        counts2 = histcounts(stationary, edges);    % Histogram for stationary data
        
        % Smooth the histograms using specified bin count for smoother distributions
        smooth_counts1 = smooth(counts1, smooth_bins); 
        smooth_counts2 = smooth(counts2, smooth_bins);
        
        % Find the peak (maximum) index for each smoothed distribution
        [~, peak1_idx] = max(smooth_counts1);  % Peak index for motion distribution
        [~, peak2_idx] = max(smooth_counts2);  % Peak index for stationary distribution
        
        % Determine threshold based on relative position of stationary and motion peaks
        if peak2_idx < peak1_idx
            % If stationary peak is before motion peak, find the intersection point
            for i = peak2_idx:length(edges)-1
                if smooth_counts1(i) > smooth_counts2(i)
                    mask_threshold = edges(i);  % Set threshold at intersection point
                    return;
                end
            end
        else
            % If no intersection, use the last edge value as the threshold
            mask_threshold = edges(length(edges));
            disp('No intersection found after the peak of the second distribution.');
        end
    end
    