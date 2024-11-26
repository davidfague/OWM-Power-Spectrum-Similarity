% Assume channel_results is a cell array of structs
% Initialize the output struct
clusters = struct();
clusters.anat = {};
clusters.num_sig_test_clusters = {};
clusters.sig_test_cluster_sizes = {};
clusters.real_test_cluster_sizes = {};

% Iterate through the cell array and concatenate fields
for chan_idx = 1:numel(channel_results)
    if isempty(channel_results{chan_idx})
        continue; % Skip empty cells
    end

    current_struct = channel_results{chan_idx};

    % Append data to the new struct fields
    clusters.anat = [clusters.anat, current_struct.anat];
    clusters.num_sig_test_clusters = [clusters.num_sig_test_clusters, current_struct.num_sig_test_clusters];
    clusters.sig_test_cluster_sizes = [clusters.sig_test_cluster_sizes, current_struct.sig_test_cluster_sizes];
    clusters.real_test_cluster_sizes = [clusters.real_test_cluster_sizes, current_struct.real_test_cluster_sizes];
end
