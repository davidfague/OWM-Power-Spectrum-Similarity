function plot_horizontal_means(data_vector)
    % Use the same index ranges as defined in add_event_lines function
    index_adjust = -5;  % Center adjustment (must match the adjustment in add_event_lines)
    
    % Define the ranges based on the event lines in add_event_lines
    index_ranges = [
        25 + index_adjust, 75 + index_adjust;  % Baseline period
        75 + index_adjust, 100 + index_adjust; % Transition from baseline to fixation
        100 + index_adjust, 150 + index_adjust; % Fixation period to end of enc1
        150 + index_adjust, 200 + index_adjust; % Enc1 to end of enc2
        200 + index_adjust, 250 + index_adjust; % Enc2 to end of enc3
        250 + index_adjust, 640 + index_adjust; % Enc3 to end of maintenance
    ];
    
    % Check if the figure is held, if not hold it
    wasHold = ishold;
    hold on;
    
    % Loop through the specified index ranges
    for i = 1:size(index_ranges, 1)
        % Get the current range
        start_idx = index_ranges(i, 1);
        end_idx = index_ranges(i, 2);
        
        % Calculate mean value for the range
        mean_value = mean(data_vector(start_idx:end_idx));
        
        % Plot a horizontal line at the mean value within the given range
        plot([start_idx, end_idx], [mean_value, mean_value], 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
    end
    
    % Restore hold state
    if ~wasHold
        hold off;
    end
end
