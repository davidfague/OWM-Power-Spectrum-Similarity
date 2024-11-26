function visualize_significant_clusters(real_p_val_matrix, real_cluster_sizes, threshold, alpha)
    figure;
    imagesc(real_p_val_matrix);
    hold on;
    colormap('gray');
    h = colorbar;
    ylabel(h, 'P-Value');
    caxis([0 alpha]);

    % Find connected components
    real_mask = real_p_val_matrix < alpha;
    connectedComponents = bwconncomp(real_mask);
    significant_clusters = real_cluster_sizes > threshold;
    
    % Plot clusters
    for i = 1:length(real_cluster_sizes)
        cluster_indices = connectedComponents.PixelIdxList{i};
        [row, col] = ind2sub(size(real_p_val_matrix), cluster_indices);
        
        if significant_clusters(i)
            plot(col, row, 'r.', 'MarkerSize', 10);
        else
            plot(col, row, 'b.', 'MarkerSize', 5);
        end
    end
    
    title('Clusters with Significant Size');
    xlabel('Maintenance Indices');
    ylabel('Encoding Indices');
    legend('Significant Cluster', 'Non-Significant Cluster', 'Location', 'Best');
    grid on;
    hold off;
end