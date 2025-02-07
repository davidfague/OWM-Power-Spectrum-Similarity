function plot_horizontal_means(data_vector, options)
    % plot_horizontal_means - Plots horizontal means with optional line style and color
    % 
    % Inputs:
    %   data_vector: Vector of data values
    %   options (optional): Structure with the following fields:
    %       - 'DashedLines' (logical): Use dashed lines if true (default: false)
    %       - 'LineColor' (1x3 vector): RGB color for the lines (default: [0.6 0.6 0.6])
    %
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

    % Set default options if not provided
    if nargin < 2 || isempty(options)
        options.DashedLines = false;
        options.LineColor = [0.6 0.6 0.6];
    end
    
    % Set line style based on options
    if options.DashedLines
        lineStyle = '--';
    else
        lineStyle = '-';
    end
    
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
        plot([start_idx, end_idx], [mean_value, mean_value], ...
            'LineWidth', 2, 'Color', options.LineColor, 'LineStyle', lineStyle);
    end
    
    % Restore hold state
    if ~wasHold
        hold off;
    end
end

