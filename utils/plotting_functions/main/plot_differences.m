function fig = plot_differences(t_test_info, plot_params)
% plot_differences plots the distribution of surrogate differences and overlays
% a vertical line at the mean of real differences. It also computes the tail
% percentage: the percentage of surrogate means that lie beyond the real 
% difference mean (to the left if the real mean is lower than the surrogate mean,
% or to the right if it is higher).
%
% Input:
%   t_test_info - structure with fields:
%       real_diff: [1×1×1000 double]
%       surrogate_diffs: [1×1×1000×1000 double]
%
% Example:
%   plot_differences(t_test_info)

    % Compute the mean of surrogate_diffs across the 3rd dimension.
    % The result is a [1 x 1 x 1 x 1000] array, so squeeze to obtain a vector.
    surrogate_means = squeeze(mean(t_test_info.surrogate_diffs, 3));  % 1x1000 vector

    % Compute the mean of real_diff across the 3rd dimension.
    real_mean = mean(t_test_info.real_diff, 3);
    real_mean = real_mean(:);    % ensure it is a vector
    real_mean = real_mean(1);    % extract the scalar value

    % Compute the overall mean of the surrogate distribution.
    surrogate_overall_mean = mean(surrogate_means);

    % Determine which tail to consider.
    if real_mean < surrogate_overall_mean
        % If real_mean is less than the surrogate overall mean, consider the left tail.
        tail_percentage = sum(surrogate_means <= real_mean) / numel(surrogate_means) * 100;
        tail_direction = 'left';
    else
        % If real_mean is greater, consider the right tail.
        tail_percentage = sum(surrogate_means >= real_mean) / numel(surrogate_means) * 100;
        tail_direction = 'right';
    end

    % Create a new figure for the plot.
    fig = figure;

    % Plot the histogram of surrogate means (normalized to form a PDF).
    h_hist = histogram(surrogate_means, 'Normalization', 'pdf');
    hold on;

    % Get current y-axis limits to set the vertical line height.
    y_limits = ylim;
    h_line = plot([real_mean, real_mean], y_limits, 'r', 'LineWidth', 2);

    % Label the axes and add a title.
    xlabel('Difference Value');
    ylabel('Probability Density');
    title_str = sprintf("Real difference on surrogate distribution\n%s p%s chan%s image%s enc%s\n %s", ...
        plot_params.type, num2str(plot_params.patient_id), ... 
        num2str(plot_params.chan_id), num2str(plot_params.image_id), ...
        num2str(plot_params.enc_id), plot_params.anat);
    title(title_str);

    % Create a legend entry for the vertical line including tail percentage info.
    legend_str = sprintf('Real Diff Mean (Tail = %.2f%%, %s)', tail_percentage, tail_direction);
    legend([h_hist, h_line], {'Surrogate Differences', legend_str}, 'Location', 'Best');

    hold off;
end
