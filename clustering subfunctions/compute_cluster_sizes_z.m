function cluster_sizes = compute_cluster_sizes_z(z_val_matrix, z_thresh)
    mask = z_val_matrix > z_thresh;
    connectedComponents = bwconncomp(mask);
    cluster_sizes = cellfun(@numel, connectedComponents.PixelIdxList);
end