%% Initialization Parameters
hellbender = false;
within_item_vs_btwn = true;
nPermutations = 1000;
use_z = true; % false uses p
alpha = 0.05;
z_thresh = 1.96;

plot_cluster_sizes_by_obs = 1:3; % choose obs to plot or set false
plot_sig_clusters_sizes_by_obs = 1:3; % choose obs to plot or set false

recompute_existing_cluster_stats = true;

%% Setup Paths
if hellbender
    addpath 'Raw Data Storage'
    addpath 'subfunctions'
    addpath 'clustering subfunctions'
    output_folder = '/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc';
    patient_IDs = [201907 201908 201903 201905 201906 201901 201910 201915];
else
    addpath('../Z_Raw Data Storage')
    addpath('../subfunctions')
    addpath('../clustering subfunctions')
    output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
    patient_IDs = [201907];
end

%% Initialize Parallel Pool
if isempty(gcp('nocreate'))
    if hellbender
        parpool('local', 64);
    else
        parpool('local', 2);
    end
end

%% Main Analysis Loop
clusters_all_patients = struct();
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    fprintf("Processing patient %d...\n", patient_ID);
    
    [label_table, all_matrix] = load_patient_data(output_folder, patient_ID);
    clusters_all_patients(patient_ID) = process_patient_data(label_table, all_matrix, nPermutations, alpha, use_z, z_thresh, recompute_existing_cluster_stats);
end

%% Plot Results Across Anatomies
plot_clusters_across_anatomies(clusters_all_patients);

%% Function Definitions

function [label_table, all_matrix] = load_patient_data(output_folder, patient_ID)
    PS_file = get_PS_file(output_folder, patient_ID);
    test_file = get_ES_file(PS_file);
    
    label_table = PS_file.label_table;
    label_table.anatomical_label = string(label_table.anatomical_label);
    rows_without_nan = get_rows_without_nan(label_table);
    label_table = label_table(rows_without_nan, :);
    
    all_matrix = squeeze(test_file.all3_ES_matrix(rows_without_nan, :, :, :));
end

function clusters_by_patient = process_patient_data(label_table, all_matrix, nPermutations, alpha, use_z, z_thresh, recompute_existing_cluster_stats)
    clusters_by_patient = struct();
    unique_channel_IDs = unique(label_table.channel_ID);
    
    for target_enc_ids = 1:size(all_matrix, 4)
        for target_image_ids = 1:size(all_matrix, 3)
            fprintf("Encoding ID: %d, Image ID: %d\n", target_enc_ids, target_image_ids);
            target_matrix = squeeze(all_matrix(:, :, target_image_ids, target_enc_ids));
            
            % Pre-slice data for parallel processing
            [pre_sliced_tables, pre_sliced_matrices] = pre_slice_data(label_table, target_matrix, unique_channel_IDs);
            
            % Compute clusters
            clusters_by_patient(target_enc_ids, target_image_ids) = compute_clusters(pre_sliced_tables, pre_sliced_matrices, ...
                unique_channel_IDs, nPermutations, alpha, use_z, z_thresh, recompute_existing_cluster_stats);
        end
    end
end

function [pre_sliced_tables, pre_sliced_matrices] = pre_slice_data(label_table, target_matrix, unique_channel_IDs)
    num_channels = length(unique_channel_IDs);
    pre_sliced_tables = cell(1, num_channels);
    pre_sliced_matrices = cell(1, num_channels);
    
    for chan_idx = 1:num_channels
        chan_id = unique_channel_IDs(chan_idx);
        chan_rows = label_table.channel_ID == chan_id;
        pre_sliced_tables{chan_idx} = label_table(chan_rows, :);
        pre_sliced_matrices{chan_idx} = target_matrix(chan_rows, :, :);
    end
end

function clusters = compute_clusters(pre_sliced_tables, pre_sliced_matrices, unique_channel_IDs, nPermutations, alpha, use_z, z_thresh, recompute)
    clusters = initialize_clusters_struct(length(unique_channel_IDs));
    parfor chan_idx = 1:length(unique_channel_IDs)
        chan_test_table = pre_sliced_tables{chan_idx};
        chan_test_matrix = squeeze(mean(pre_sliced_matrices{chan_idx}, 1));
        control_matrices = generate_control_matrices(chan_test_matrix, nPermutations);
        
        % Perform clustering analysis
        [num_clusters, sig_sizes, real_sizes, ~, threshold] = ...
            get_clusters(control_matrices, chan_test_matrix, nPermutations, alpha, use_z, z_thresh);
        
        % Store results
        anat = unique(chan_test_table.anatomical_label);
        clusters = store_cluster_results(clusters, chan_idx, anat, num_clusters, sig_sizes, real_sizes);
    end
end

function plot_clusters_across_anatomies(clusters_all_patients)
    anat_labels = [];
    cluster_sizes = [];
    
    % Aggregate data
    patient_IDs = fieldnames(clusters_all_patients);
    for i = 1:length(patient_IDs)
        patient_clusters = clusters_all_patients.(patient_IDs{i});
        anat_labels = [anat_labels; patient_clusters.anat];
        cluster_sizes = [cluster_sizes; cell2mat(patient_clusters.real_test_cluster_sizes)];
    end
    
    % Plot aggregated data
    unique_labels = unique(anat_labels);
    mean_cluster_sizes = arrayfun(@(x) mean(cluster_sizes(anat_labels == x)), unique_labels);
    
    figure;
    bar(categorical(unique_labels), mean_cluster_sizes);
    xlabel('Anatomical Label');
    ylabel('Mean Cluster Size');
    title('Cluster Sizes Across Anatomies');
end
