function concatenatedVector = concatenate_cluster_sizes(clusterSizesNull)
    concatenatedVector = nan(sum(cellfun(@numel, clusterSizesNull)), 1);
    for i = 1:length(clusterSizesNull)
        concatenatedVector = [concatenatedVector; clusterSizesNull{i}(:)]; % Ensure column vector
    end
end