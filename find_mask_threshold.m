
function [mask_threshold, smooth_counts1, smooth_counts2] = find_mask_threshold(motion, stationary, edges)
    % Inputs:
    %   motion:
    %   stationary:
    %   edges: bin edges for the histogram
    
    % Compute histograms
    counts1 = histcounts(motion, edges);
    counts2 = histcounts(stationary, edges);
    
    % Smooth the histograms
    smooth_counts1 = smooth(counts1, 100); 
    smooth_counts2 = smooth(counts2, 100);
    
    % Find the peaks of the distributions
    [~, peak1_idx] = max(smooth_counts1);
    [~, peak2_idx] = max(smooth_counts2);
    
    % If the peak of the stationary distribution comes before of the
    % motion peak, we set a threshold
    if peak2_idx < peak1_idx
        % Find the first intersection point after the peak of the second distribution
        for i = peak2_idx:length(edges)-1
            if smooth_counts1(i) > smooth_counts2(i)
                mask_threshold = edges(i);
                return;
            end
        end
    else
        % Else we take all the data
        mask_threshold = edges(length(edges));
        disp('No intersection found after the peak of the second distribution.');
    end
end