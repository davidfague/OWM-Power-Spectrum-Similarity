function fig = plot_cluster_sum_distribution_copy(all_surrogate_sums, all_real_cluster_sums, significance_threshold, final_p_this_obs)
    sum_threshold = prctile(abs(all_surrogate_sums), (1-significance_threshold)*100);
 % Generate kernel density estimates for each dataset
    [density_surrogate, x_surrogate] = ksdensity(abs(all_surrogate_sums));
    [density_real, x_real] = ksdensity(all_real_cluster_sums);
    
    % Plot surrogate distributions as smooth lines
    fig = figure;
    plot(x_surrogate, density_surrogate, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(x_real, density_real, 'k-', 'LineWidth', 1.5);

    
    % Add thresholds as vertical lines
    xline(sum_threshold, 'r--', 'LineWidth', 2, 'Label', sprintf('p < %s', num2str(significance_threshold)));
    xline(max(all_real_cluster_sums), 'bl--', 'LineWidth',2, 'Label', sprintf('p = %s', num2str(final_p_this_obs)))

    % Add labels, title, and legend
    title('Time Cluster Sizes by Sum(t-stat) (PDF)');
    xlabel('Cluster Sum(t-stat)');
    ylabel('Probability Density');
    legend({'Null Nperm=1000', 'Real', ...
        sprintf('Significance Threshold (p < %s)', num2str(significance_threshold)), ...
        sprintf('Real Cluster Max (p < %s)', num2str(final_p_this_obs))}, ...
        'Location', 'best');
    % hold off;
end