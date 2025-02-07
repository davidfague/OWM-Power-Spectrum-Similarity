
function visualize_significant_clusters(real_p_val_matrix, real_cluster_sizes, threshold, alpha)
    figure;
    imagesc(real_p_val_matrix);
    hold on;
    colormap('gray');
    h = colorbar;
    ylabel(h, 'P-Value');
    clim([0 alpha]);

    % Find connected components
    real_mask = real_p_val_matrix < alpha;
    connectedComponents = bwconncomp(real_mask);
    significant_clusters = real_cluster_sizes > threshold;

    % Initialize flags for legend
    plotted_significant = false;
    plotted_non_significant = false;

    % Plot clusters
    for i = 1:length(real_cluster_sizes)
        cluster_indices = connectedComponents.PixelIdxList{i};
        [row, col] = ind2sub(size(real_p_val_matrix), cluster_indices);

        if significant_clusters(i)
            plot(col, row, 'r.', 'MarkerSize', 10);
            plotted_significant = true;
        else
            plot(col, row, 'b.', 'MarkerSize', 5);
            plotted_non_significant = true;
        end
    end

    title('Clusters with Significant Size');
    xlabel('Maintenance Indices');
    ylabel('Encoding Indices');

    % Adjust legend based on plotted clusters
    if plotted_significant && plotted_non_significant
        legend('Significant Cluster', 'Non-Significant Cluster', 'Location', 'Best');
    elseif plotted_significant
        legend('Significant Cluster', 'Location', 'Best');
    elseif plotted_non_significant
        legend('Non-Significant Cluster', 'Location', 'Best');
    end

    grid on;
    hold off;
end