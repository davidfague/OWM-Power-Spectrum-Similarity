function plot_cluster_sizes_by_uniq_anat(clusters_by_uniq_anat)
    unique_anat = clusters_by_uniq_anat.anat;
    n_anat = length(unique_anat);
    cluster_sizes = clusters_by_uniq_anat.real_test_cluster_sizes;

    all_cluster_sizes = [];
    group_labels = [];
    for anat_idx = 1:n_anat
        all_cluster_sizes = [all_cluster_sizes; cluster_sizes{anat_idx}];
        group_labels = [group_labels; repmat(unique_anat(anat_idx), length(cluster_sizes{anat_idx}), 1)];
    end

    figure;
    boxplot(all_cluster_sizes, group_labels, 'LabelOrientation', 'inline');
    xlabel('Anatomical Labels');
    ylabel('Cluster Size');
    title('Cluster Sizes Grouped by Anatomical Regions');
    set(gca, 'XTickLabelRotation', 45);
    grid on;
end