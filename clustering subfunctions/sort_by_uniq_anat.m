function clusters_by_uniq_anat = sort_by_uniq_anat(clusters)
    unique_anat = unique(clusters.anat);
    n_anat = length(unique_anat);

    clusters_by_uniq_anat = struct();
    clusters_by_uniq_anat.anat = unique_anat;
    clusters_by_uniq_anat.num_sig_test_clusters = cell(n_anat, 1);
    clusters_by_uniq_anat.sig_test_cluster_sizes = cell(n_anat, 1);
    clusters_by_uniq_anat.real_test_cluster_sizes = cell(n_anat, 1);

    for anat_idx = 1:n_anat
        anat_label = unique_anat{anat_idx};
        matching_indices = strcmp(clusters.anat, anat_label);

        clusters_by_uniq_anat.num_sig_test_clusters{anat_idx} = ...
            vertcat(clusters.num_sig_test_clusters(matching_indices));
        clusters_by_uniq_anat.sig_test_cluster_sizes{anat_idx} = ...
            vertcat(clusters.sig_test_cluster_sizes{matching_indices});
        clusters_by_uniq_anat.real_test_cluster_sizes{anat_idx} = ...
            vertcat(clusters.real_test_cluster_sizes{matching_indices});
    end
end