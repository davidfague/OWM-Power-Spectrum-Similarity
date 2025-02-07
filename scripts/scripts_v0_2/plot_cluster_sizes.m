%% Collect and Aggregate Cluster Sizes
all_cluster_sizes = struct();
all_cluster_sizes.anat = {};
all_cluster_sizes.real_test_cluster_sizes = {};

for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    PS_file = get_PS_file(output_folder, patient_ID);
    test_file = get_ES_file(PS_file);
    label_table = PS_file.label_table;
    label_table.anatomical_label = string(label_table.anatomical_label);
    rows_without_nan = get_rows_without_nan(label_table);
    label_table = label_table(rows_without_nan, :);

    all_matrix = squeeze(test_file.all3_ES_matrix(:, :, :, :));
    all_matrix = all_matrix(rows_without_nan, :, :, :);

    for target_enc_ids = 1:3 % Loop through encoding IDs
        target_matrix = squeeze(all_matrix(:, :, :, target_enc_ids));

        for target_image_ids = 1:9 % Loop through image IDs
            ctrl_file = get_ctrl_file(PS_file, patient_ID, target_enc_ids, target_image_ids);
            unique_channel_IDs = ctrl_file.unique_channel_IDs;

            clusters_save_file = strrep(ctrl_file.Properties.Source, 'all_channels_data.mat', 'all_channel_cluster_results.mat');
            if exist(clusters_save_file, 'file')
                data = load(clusters_save_file, "channel_results");
                channel_results = data.channel_results;
                clear data

                % Process and collect clusters
                clusters = fix_structure(channel_results);
                all_cluster_sizes.anat = [all_cluster_sizes.anat; clusters.anat];
                all_cluster_sizes.real_test_cluster_sizes = [all_cluster_sizes.real_test_cluster_sizes; clusters.real_test_cluster_sizes];
            end
        end
    end
end

%% Aggregate Cluster Sizes by Unique Anatomical Labels
unique_anats = unique(all_cluster_sizes.anat);
cluster_size_by_anat = cell(size(unique_anats));

for i = 1:length(unique_anats)
    anat = unique_anats{i};
    mask = strcmp(all_cluster_sizes.anat, anat);
    cluster_sizes = vertcat(all_cluster_sizes.real_test_cluster_sizes{mask});
    cluster_size_by_anat{i} = cluster_sizes;
end

%% Plot Cluster Sizes
figure;
hold on;

for i = 1:length(unique_anats)
    anat = unique_anats{i};
    sizes = cluster_size_by_anat{i};
    
    % Use a jittered x position to better visualize overlapping points
    x = i + 0.1 * randn(size(sizes));
    scatter(x, sizes, 'filled');
end

xticks(1:length(unique_anats));
xticklabels(unique_anats);
xtickangle(45);
xlabel('Anatomical Label');
ylabel('Cluster Size');
title('Cluster Sizes by Anatomical Label');
grid on;
hold off;


