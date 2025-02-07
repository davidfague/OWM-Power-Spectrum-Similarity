function cluster_analysis_and_visualization_abs(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values)
    % Parameters
    alpha = 0.025; % Significance level
    significance_threshold = 0.025; % For cluster significance

    % Step 1: Threshold the real data for positive and negative clusters
    positive_mask = (real_p_values < alpha) & (real_t_values > 0);
    negative_mask = (real_p_values < alpha) & (real_t_values < 0);

    % % Plot thresholds
    % figure; hold on;
    % subplot(1, 2, 1);
    % imagesc(positive_mask);
    % title('Positive Clusters (p<0.05, t>0)');
    % xlabel('yTimes');
    % ylabel('xTimes');
    % colorbar;
    % 
    % subplot(1, 2, 2);
    % imagesc(negative_mask);
    % title('Negative Clusters (p<0.05, t<0)');
    % xlabel('yTimes');
    % ylabel('xTimes');
    % colorbar;

    % Step 2: Identify connected components (clusters)
    positive_clusters = bwconncomp(positive_mask, 8); % 8-connectivity
    negative_clusters = bwconncomp(negative_mask, 8);

    % Step 3: Calculate cluster sums for real data
    real_positive_sums = calculate_cluster_sums(real_t_values, positive_clusters);
    real_negative_sums = calculate_cluster_sums(real_t_values, negative_clusters);

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
    positive_p_values = assess_significance(real_positive_sums, surrogate_positive_sums, 'positive');
    negative_p_values = assess_significance(real_negative_sums, surrogate_negative_sums, 'negative');

    % Determine significance for clusters
    significant_positive = positive_p_values < significance_threshold;
    significant_negative = negative_p_values < significance_threshold;

    % compute the threshold sum for plotting
    positive_sum_threshold = prctile(surrogate_positive_sums, (1-significance_threshold)*100);
    negative_sum_threshold = prctile(surrogate_negative_sums, significance_threshold*100);

    % Step 6: Visualization
    plot_cluster_sum_distribution(surrogate_positive_sums, surrogate_negative_sums, real_negative_sums, real_positive_sums, significance_threshold)
    %%
    % Visualize detected clusters
    figure; hold on;
    subplot(1,2,1)
    visualize_clusters(real_p_values, real_t_values, positive_clusters, significant_positive, positive_sum_threshold, alpha, 'Positive');
    subplot(1,2,2)
    visualize_clusters(real_p_values, real_t_values, negative_clusters, significant_negative, negative_sum_threshold, alpha, 'Negative');
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

function plot_cluster_sum_distribution(surrogate_positive_sums, surrogate_negative_sums, real_negative_sums, real_positive_sums, significance_threshold)
    positive_sum_threshold = prctile(surrogate_positive_sums, (1-significance_threshold)*100);
    negative_sum_threshold = prctile(surrogate_negative_sums, significance_threshold*100);
 % Generate kernel density estimates for each dataset
    [density_surrogate_positive, x_surrogate_positive] = ksdensity(surrogate_positive_sums);
    [density_surrogate_negative, x_surrogate_negative] = ksdensity(surrogate_negative_sums);
    % [density_real_positive, x_real_positive] = ksdensity(real_positive_sums);
    % [density_real_negative, x_real_negative] = ksdensity(real_negative_sums);
    
    % Plot surrogate distributions as smooth lines
    figure; hold on;
    h1 = plot(x_surrogate_positive, density_surrogate_positive, 'b-', 'LineWidth', 1.5);
    h2 = plot(x_surrogate_negative, density_surrogate_negative, 'r-', 'LineWidth', 1.5);
    % plot(x_real_positive, density_real_positive, 'cyan--', 'LineWidth', 1.5);
    % plot(x_real_negative, density_real_negative, 'o--', 'LineWidth', 1.5);

    for i=1:length(real_negative_sums)
        xline(real_negative_sums(i), 'Color','magenta')%, 'Label','Real Negative Clusters')
    end

    for i=1:length(real_positive_sums)
        xline(real_positive_sums(i), 'Color','cyan')%, 'Label','Real Positive Clusters')
    end

    
    % Add thresholds as vertical lines
    xline(positive_sum_threshold, 'b--', 'LineWidth', 2, 'Label', 'Positive Threshold');
    xline(negative_sum_threshold, 'r--', 'LineWidth', 2, 'Label', 'Negative Threshold','LabelHorizontalAlignment', 'left');
    
    % Create dummy handles for legend
    h3 = plot(NaN, NaN, 'c-', 'LineWidth', 1.5); % Real Positive Clusters
    h4 = plot(NaN, NaN, 'm-', 'LineWidth', 1.5); % Real Negative Clusters
    h5 = plot(NaN, NaN, 'b--', 'LineWidth', 2); % Positive Threshold
    h6 = plot(NaN, NaN, 'r--', 'LineWidth', 2); % Negative Threshold

    % Add labels, title, and legend
    title('Positive & Negative Clusters - Surrogate Distribution');
    xlabel('Cluster Sum');
    ylabel('Density');
    legend([h1, h2, h3, h4, h5, h6], ...
        {'Surrogate Positive', 'Surrogate Negative', ...
         'Real Positive Clusters', 'Real Negative Clusters',...
        sprintf('Positive Significance Threshold (p < %s)', num2str(significance_threshold)), ...
        sprintf('Negative Significance Threshold (p < %s)', num2str(significance_threshold))}, ...
        'Location', 'best');
    hold off;
    end

function visualize_clusters(p_values, t_values, clusters, significance, threshold, threshold_alpha, cluster_type)
    % Visualize clusters and label them as significant or not
    
    % Initialize cluster_mask with the same size as t_values
    cluster_mask = false(size(p_values));
    
    % Populate cluster_mask based on clusters
    for i = 1:clusters.NumObjects
        cluster_mask(clusters.PixelIdxList{i}) = true; % Mark cluster indices
    end

    p_values(~cluster_mask) = 1; % effectively remove p values that are not in the clusters of this type.

    % Plot the grayscale background
    % figure;
    imagesc(p_values);
    colormap('gray');
    clim([0, threshold_alpha]); % Ensure consistent scaling
    colorbar;
    title(sprintf('%s Clusters (Threshold = %.2f)', cluster_type, threshold));
    xlabel('Maintenance Time Bins (10 ms)');
    ylabel('Enc Time Bins (10 ms)');
    hold on;

    % Overlay clusters with red or blue dots based on significance
    plotted_significant = false;
    plotted_non_significant = false;
    for i = 1:clusters.NumObjects
        [rows, cols] = ind2sub(size(p_values), clusters.PixelIdxList{i});

        % Calculate the sum of t_values for the current cluster
        sum_t = sum(t_values(clusters.PixelIdxList{i}));

        % Determine the cluster's position for labeling (mean of cluster coordinates)
        label_row = mean(rows);
        label_col = mean(cols);

        if significance(i)
            plot(cols, rows, 'r.', 'MarkerSize', 10); % Red for significant
            plotted_significant = true;

            % Add text label for the sum t value
            text(label_col, label_row, sprintf('%.2f', sum_t), 'Color', 'white', ...
                'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

        else
            plot(cols, rows, 'b.', 'MarkerSize', 5); % Blue for non-significant
            plotted_non_significant = true;

            % Add text label for the sum t value
            text(label_col, label_row, sprintf('%.2f', sum_t), 'Color', 'blue', ...
                'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end

    end

    % Adjust legend
    if plotted_significant && plotted_non_significant
        legend('Significant Cluster', 'Non-Significant Cluster', 'Location', 'Best');
    elseif plotted_significant
        legend('Significant Cluster', 'Location', 'Best');
    elseif plotted_non_significant
        legend('Non-Significant Cluster', 'Location', 'Best');
    end

    % hold off;
end
