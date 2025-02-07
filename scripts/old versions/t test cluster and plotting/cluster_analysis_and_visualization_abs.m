function [final_p_this_obs] = cluster_analysis_and_visualization_abs(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values)
    % Parameters
    alpha = 0.025; % Significance level For separating clusters
    significance_threshold = 0.05; % For cluster T sum significance

    % Step 1: Threshold the real data for positive and negative clusters
    positive_mask = (real_p_values < alpha) & (real_t_values > 0);
    negative_mask = (real_p_values < alpha) & (real_t_values < 0);

    % % Plot thresholds
    figure; hold on;
    subplot(1, 2, 1);
    imagesc(positive_mask);
    title('Positive Clusters (p<0.05, t>0)');
    xlabel('yTimes');
    ylabel('xTimes');
    colormap('jet')
    colorbar;

    subplot(1, 2, 2);
    imagesc(negative_mask);
    title('Negative Clusters (p<0.05, t<0)');
    xlabel('yTimes');
    ylabel('xTimes');
    colormap('jet')
    colorbar;

    figure;
    subplot(1,2,1)
    imagesc(real_t_values)
    colormap('jet')
    colorbar;
    clim([-30, 30])
    subplot(1,2,2)
    imagesc(surrogate_t_values(:,:,1))
    colormap('jet')
    colorbar;
    clim([-30, 30])

    % Step 2: Identify connected components (clusters)
    positive_clusters = bwconncomp(positive_mask, 8); % 8-connectivity
    negative_clusters = bwconncomp(negative_mask, 8);
    all_real_clusters = [positive_clusters, negative_clusters];

    % Step 3: Calculate cluster sums for real data
    real_positive_sums = calculate_cluster_sums(real_t_values, positive_clusters);
    real_negative_sums = calculate_cluster_sums(real_t_values, negative_clusters);
    all_real_cluster_sums = [real_positive_sums, abs(real_negative_sums)];

    % Step 4: Compute surrogate cluster sums for the largest cluster
    num_surrogates = size(surrogate_t_values, 3);
    surrogate_positive_sums = zeros(num_surrogates, 1);
    surrogate_negative_sums = zeros(num_surrogates, 1);
    all_surrogate_sums = zeros(num_surrogates, 1);

    for i = 1:num_surrogates
        % Threshold surrogate data
        surrogate_positive_mask = (surrogate_p_values(:, :, i) < alpha) & (surrogate_t_values(:, :, i) > 0);
        surrogate_negative_mask = (surrogate_p_values(:, :, i) < alpha) & (surrogate_t_values(:, :, i) < 0);

        % Find clusters
        surrogate_positive_clusters = bwconncomp(surrogate_positive_mask, 8);
        surrogate_negative_clusters = bwconncomp(surrogate_negative_mask, 8);

        % Calculate largest cluster sum
        surrogate_positive_sums(i) = max(calculate_cluster_sums(surrogate_t_values(:, :, i), surrogate_positive_clusters), [], 'omitnan');
        surrogate_negative_sums(i) = min(calculate_cluster_sums(surrogate_t_values(:, :, i), surrogate_negative_clusters), [], 'omitnan');
        if abs(surrogate_positive_sums(i)) > abs(surrogate_negative_sums(i))
            all_surrogate_sums(i) = surrogate_positive_sums(i);
        else
            all_surrogate_sums(i) = surrogate_negative_sums(i);
        end
    end

    % Step 5: Assess significance of real clusters
    % positive_p_values = assess_significance(real_positive_sums, surrogate_positive_sums, 'positive');
    % negative_p_values = assess_significance(real_negative_sums, surrogate_negative_sums, 'negative');

    cluster_p_values = assess_significance(all_real_cluster_sums, all_surrogate_sums, 'positive');

    % Determine significance for clusters
    % significant_positive = positive_p_values < significance_threshold;
    % significant_negative = negative_p_values < significance_threshold;
    significant_clusters = cluster_p_values < significance_threshold;

    % compute the threshold sum for plotting
    % positive_sum_threshold = prctile(surrogate_positive_sums, (1-significance_threshold)*100);
    % negative_sum_threshold = prctile(surrogate_negative_sums, significance_threshold*100);
    sum_threshold = prctile(abs(all_surrogate_sums), (1-significance_threshold)*100);

    final_p_this_obs = mean(abs(all_surrogate_sums)>max(all_real_cluster_sums));

    % Step 6: Visualization
    % plot_cluster_sum_distribution(surrogate_positive_sums, surrogate_negative_sums, real_negative_sums, real_positive_sums, significance_threshold)
    plot_cluster_sum_distribution(all_surrogate_sums, all_real_cluster_sums, significance_threshold)
    %%
    % Visualize detected clusters
    % figure; hold on;
    % subplot(1,2,1)
    % visualize_clusters(real_p_values, real_t_values, positive_clusters, significant_positive, positive_sum_threshold, alpha, 'Positive');
    % subplot(1,2,2)
    % visualize_clusters(real_p_values, real_t_values, negative_clusters, significant_negative, negative_sum_threshold, alpha, 'Negative');
    visualize_clusters(real_p_values, real_t_values, all_real_clusters, significant_clusters, sum_threshold, alpha);
end

function cluster_sums = calculate_cluster_sums(t_values, clusters)
    % Calculate sum of t-values for each cluster
    cluster_sums = zeros(1, clusters.NumObjects);
    for i = 1:clusters.NumObjects
        cluster_sums(i) = sum(t_values(clusters.PixelIdxList{i}), 'omitnan');
    end
end

function p_values = assess_significance(real_sums, surrogate_sums, cluster_type)
    % Compute p-values by comparing real cluster sums to surrogate distribution
    p_values = zeros(size(real_sums));
    if strcmp(cluster_type, 'positive')
        % Positive clusters: Compare to max surrogate sums
        for i = 1:length(real_sums)
            p_values(i) = mean(surrogate_sums >= real_sums(i));
        end
    elseif strcmp(cluster_type, 'negative')
        % Negative clusters: Compare to min surrogate sums
        for i = 1:length(real_sums)
            p_values(i) = mean(surrogate_sums <= real_sums(i));
        end
    else
        error('Invalid cluster type. Use "positive" or "negative".');
    end
end

function plot_cluster_sum_distribution(all_surrogate_sums, all_real_cluster_sums, significance_threshold)
    sum_threshold = prctile(abs(all_surrogate_sums), (1-significance_threshold)*100);
 % Generate kernel density estimates for each dataset
    [density_surrogate, x_surrogate] = ksdensity(abs(all_surrogate_sums));
    [density_real, x_real] = ksdensity(all_real_cluster_sums);
    
    % Plot surrogate distributions as smooth lines
    figure;
    plot(x_surrogate, density_surrogate, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(x_real, density_real, 'k-', 'LineWidth', 1.5);

    
    % Add thresholds as vertical lines
    xline(sum_threshold, 'r--', 'LineWidth', 2, 'Label', sprintf('p < %s', num2str(significance_threshold)));

    % Add labels, title, and legend
    title('Positive & Negative Clusters - Surrogate Distribution');
    xlabel('Cluster Sum');
    ylabel('Density');
    legend({'Null Nperm=1000', 'Real', ...
        sprintf('Significance Threshold (p < %s)', num2str(significance_threshold))}, ...
        'Location', 'best');
    % hold off;
end

function visualize_clusters(p_values, t_values, clusters, significance, threshold, threshold_alpha)
    % Visualize clusters and label them as significant or not
    
    % Initialize cluster_mask with the same size as t_values
    % cluster_mask = false(size(p_values));
    
    % % Populate cluster_mask based on clusters
    % for i = 1:clusters.NumObjects
    %     cluster_mask(clusters.PixelIdxList{i}) = true; % Mark cluster indices
    % end
    % 
    % p_values(~cluster_mask) = 1; % effectively remove p values that are not in the clusters of this type.

    % Plot the grayscale background
    figure;
    imagesc(t_values);
    % colormap('gray');
    % clim([0, threshold_alpha]); % Ensure consistent scaling
    colormap('jet')
    colorbar;
    title(sprintf('Clusters (Threshold = %.2f)', threshold));
    xlabel('Maintenance Time Bins (10 ms)');
    ylabel('Enc Time Bins (10 ms)');
    hold on;

    % Overlay clusters with red or blue dots based on significance
    % Initialize flags for plotting legend
    plotted_significant_pos = false;
    plotted_significant_neg = false;
    plotted_non_significant_pos = false;
    plotted_non_significant_neg = false;
    for cluster_type=1:2 % t=pos then t=neg
        for i = 1:clusters(cluster_type).NumObjects
            [rows, cols] = ind2sub(size(p_values), clusters(cluster_type).PixelIdxList{i});
    
            % Calculate the sum of t_values for the current cluster
            sum_t = sum(abs(t_values(clusters(cluster_type).PixelIdxList{i})));
    
            % Determine the cluster's position for labeling (mean of cluster coordinates)
            label_row = mean(rows);
            label_col = mean(cols);
    
            if cluster_type == 2
                j=i +clusters(1).NumObjects; % shift to match
            else
                j=i;
            end

            if significance(j)
                if cluster_type == 1 % pos (yellow normally)
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
                if cluster_type == 1 % pos (yellow normally)
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