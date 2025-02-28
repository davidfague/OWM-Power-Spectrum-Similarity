function cluster_visualization(t_test_info, clusters_info, params)
    % Parameters
    % alpha = 0.025; % significance level For separating clusters
    % significance_threshold = 0.025; % for cluster T sum significance

    % % Plot thresholds
    fig = figure; hold on;
    subplot(1, 2, 1);
    imagesc(clusters_info.real_positive_mask);
    title('Positive Clusters (p<0.05, t>0)');
    xlabel('yTimes');
    ylabel('xTimes');
    colormap('jet')
    colorbar;

    subplot(1, 2, 2);
    imagesc(clusters_info.real_negative_mask);
    title('Negative Clusters (p<0.05, t<0)');
    xlabel('yTimes');
    ylabel('xTimes');
    colormap('jet')
    colorbar;
    saveas(fig, sprintf("%s/Raw clusters.png", params.WI_BI_folder_to_save_in))

    fig = figure;
    subplot(1,3,1)
    imagesc(t_test_info.real_t_values)
    colormap('jet')
    colorbar;
    clim([-30, 30])
    title('Real')
    subplot(1,3,2)
    imagesc(t_test_info.surrogate_t_values(:,:,1))
    colormap('jet')
    colorbar;
    clim([-30, 30])
    title('1 Null')
    subplot(1,3,3)
    imagesc(mean(t_test_info.surrogate_t_values,3))
    colormap('jet')
    colorbar;
    clim([-30, 30])
    title('avg Null')
    saveas(fig, sprintf("%s/Raw t values.fig", params.WI_BI_folder_to_save_in))

    fig = plot_cluster_sum_distribution(clusters_info.surrogate_positive_sums, ...
        clusters_info.real_positive_cluster_sums, ...
        clusters_info.significance_threshold, ...
        clusters_info.p_largest_positive_cluster, ...
        "Positive");
    if ~ isempty(fig)
        saveas(fig, sprintf("%s/Positive Cluster Tsum distribution.fig", params.WI_BI_folder_to_save_in))
    end
    fig = plot_cluster_sum_distribution(clusters_info.surrogate_negative_sums, ...
        clusters_info.real_negative_cluster_sums, ...
        clusters_info.significance_threshold, ...
        clusters_info.p_largest_negative_cluster, ...
        'Negative');
    if ~isempty(fig)
        saveas(fig, sprintf("%s/Negative Cluster Tsum distribution.fig", params.WI_BI_folder_to_save_in))
    end
    %%
    % Visualize detected clusters
    % figure; hold on;
    % subplot(1,2,1)
    % visualize_clusters(real_p_values, real_t_values, positive_clusters, significant_positive, positive_sum_threshold, alpha, 'Positive');
    % subplot(1,2,2)
    % visualize_clusters(real_p_values, real_t_values, negative_clusters, significant_negative, negative_sum_threshold, alpha, 'Negative');
    fig = visualize_clusters(t_test_info, clusters_info, 'real_diff');
    saveas(fig, sprintf("%s/WI_BI_cluster_visualizastion.png", params.WI_BI_folder_to_save_in))
end

% function fig = plot_cluster_sum_distribution(all_surrogate_sums, all_real_cluster_sums, significance_threshold, final_p_this_obs)
%     sum_threshold = prctile(abs(all_surrogate_sums), (1-significance_threshold)*100);
%  % Generate kernel density estimates for each dataset
%     [density_surrogate, x_surrogate] = ksdensity(abs(all_surrogate_sums));
%     [density_real, x_real] = ksdensity(all_real_cluster_sums);
% 
%     % Plot surrogate distributions as smooth lines
%     fig = figure;
%     plot(x_surrogate, density_surrogate, 'b-', 'LineWidth', 1.5);
%     hold on;
%     plot(x_real, density_real, 'k-', 'LineWidth', 1.5);
% 
% 
%     % Add thresholds as vertical lines
%     xline(sum_threshold, 'r--', 'LineWidth', 2, 'Label', sprintf('p < %s', num2str(significance_threshold)));
%     xline(max(all_real_cluster_sums), 'bl--', 'LineWidth',2, 'Label', sprintf('p = %s', num2str(final_p_this_obs)))
% 
%     % Add labels, title, and legend
%     title('Positive & Negative Clusters - Surrogate Distribution');
%     xlabel('Cluster Sum');
%     ylabel('Density');
%     legend({'Null Nperm=1000', 'Real', ...
%         sprintf('Significance Threshold (p < %s)', num2str(significance_threshold)), ...
%         sprintf('Real Cluster Max (p < %s)', num2str(final_p_this_obs))}, ...
%         'Location', 'best');
%     % hold off;
% end

function fig = visualize_clusters(t_test_info, clusters_info, var_to_underlay)
    values_to_underlay = t_test_info.(var_to_underlay);
    dims = ndims(values_to_underlay);

    % If the matrix has more than 2 dimensions, take the mean across all
    % dimensions beyond the first and second
    if dims > 2
        values_to_underlay = mean(values_to_underlay, dims:-1:3);
    end
    
    % Plot the grayscale background
    fig = figure;
    imagesc(values_to_underlay);
    % colormap('gray');
    % clim([0, threshold_alpha]); % Ensure consistent scaling
    colormap('jet')
    c = colorbar;
    c.Label.String = var_to_underlay;
    title(sprintf('Clusters (P<%.3f, Positive Threshold = %.f, Negative Threshold = %.f)', ...
        clusters_info.significance_threshold, clusters_info.positive_sum_threshold, clusters_info.negative_sum_threshold));
    xticklabels(cellstr(string(str2double(xticklabels) * 10)));
    yticklabels(cellstr(string(str2double(yticklabels) * 10)));
    xlabel('Maintenance Time (ms)');
    ylabel('Encoding Time (ms)');
    hold on;

    % Overlay clusters with red or blue dots based on significance
    % Initialize flags for plotting legend
    plotted_significant_pos = false;
    plotted_significant_neg = false;
    plotted_non_significant_pos = false;
    plotted_non_significant_neg = false;

    %% NOTE TO UPDATE: enable the use of diffs, p's, t's. use new clusters_info and t_test_info
    cluster_types = {'positive', 'negative'};
    for cluster_type_idx=1:2 % t=pos then t=neg
        cluster_type = cluster_types{cluster_type_idx};
        clusters = clusters_info.(sprintf('real_%s_clusters',cluster_type));
        for i = 1:clusters.NumObjects
            [rows, cols] = ind2sub(size(values_to_underlay), clusters.PixelIdxList{i});
    
            % Calculate the sum of t_values for the current cluster
            %sum_t = sum(abs(t_test_info.real_t_values(clusters.PixelIdxList{i})));
            sum_t = clusters_info.(['real_', cluster_type,'_cluster_sums'])(i); % already calculated

            % Determine the cluster's position for labeling (mean of cluster coordinates)
            label_row = mean(rows);
            label_col = mean(cols);

            if clusters_info.(['significant_', cluster_type, '_clusters'])(i)
                if strcmp(cluster_type, 'positive') % pos (yellow normally)
                    plot(cols, rows, 'r.', 'MarkerSize', 10); % Red for significant positive
                    plotted_significant_pos = true;
                else % neg (blue normally)
                    plot(cols, rows, 'm.', 'MarkerSize', 10); % Magenta for significant negative
                    plotted_significant_neg = true;
                end
    
                % Add text label for the sum t value
                text(label_col, label_row, sprintf('%.2f', sum_t), 'Color', 'white', ...
                    'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
            else
                if strcmp(cluster_type, 'positive') % pos (yellow normally)
                    plot(cols, rows, 'y.', 'MarkerSize', 5); % Yellow for non-significant positive
                    plotted_non_significant_pos = true;
                else % neg (blue normally)
                    plot(cols, rows, 'b.', 'MarkerSize', 5); % Blue for non-significant negative
                    plotted_non_significant_neg = true;
                end
    
                % Add text label for the sum t value
                text(label_col, label_row, sprintf('%.2f', sum_t), 'Color', 'black', ...
                    'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
        end
    end

    % Create the legend plot objects dynamically based on the plotted flags
    hold on;
    hObjects = {}; % Cell array to store plot handles for the legend
    % Create legend entries dynamically based on which clusters were plotted
    legend_entries = {};
    
    % Combine the plotting and legend entries assignment in one block
    if plotted_significant_pos
        hObjects{end + 1} = plot(nan, nan, 'r.', 'MarkerSize', 10); % Red for significant positive
        legend_entries{end + 1} = 'Significant Positive Cluster (Red)';
    end
    if plotted_significant_neg
        hObjects{end + 1} = plot(nan, nan, 'm.', 'MarkerSize', 10); % Magenta for significant negative
        legend_entries{end + 1} = 'Significant Negative Cluster (Magenta)';
    end
    if plotted_non_significant_pos
        hObjects{end + 1} = plot(nan, nan, 'y.', 'MarkerSize', 5); % Yellow for non-significant positive
        legend_entries{end + 1} = 'Non-Significant Positive Cluster (Yellow)';
    end
    if plotted_non_significant_neg
        hObjects{end + 1} = plot(nan, nan, 'b.', 'MarkerSize', 5); % Blue for non-significant negative
        legend_entries{end + 1} = 'Non-Significant Negative Cluster (Blue)';
    end

    % Create the legend if any clusters were plotted
    if ~isempty(legend_entries)
        legend([hObjects{:}], legend_entries{:}, 'Location', 'Best');
    end
    hold off
end