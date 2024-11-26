function cluster_sizes = compute_cluster_sizes(p_val_matrix, alpha)
    mask = p_val_matrix < alpha;
    connectedComponents = bwconncomp(mask);
    cluster_sizes = cellfun(@numel, connectedComponents.PixelIdxList);
end
