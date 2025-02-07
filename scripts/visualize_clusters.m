function fig = visualize_clusters(p_values, t_values, clusters, significance, threshold, threshold_alpha, final_p_this_obs)
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
    fig = figure;
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