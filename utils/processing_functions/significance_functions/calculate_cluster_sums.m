function cluster_sums = calculate_cluster_sums(t_values, clusters)
    % Calculate sum of t-values for each cluster
    cluster_sums = zeros(1, clusters.NumObjects);
    for i = 1:clusters.NumObjects
        cluster_sums(i) = sum(t_values(clusters.PixelIdxList{i}), 'omitnan');
    end
end