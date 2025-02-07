%% Initialization Parameters
hellbender = false;
within_item_vs_btwn = false;
nPermutations = 1000;
use_z = true; % false uses p
alpha = 0.05;
z_thresh = 1.96;

plot_cluster_sizes_by_obs = 1:3; % choose obs to plot or set false
plot_sig_clusters_sizes_by_obs = 1:3; % choose obs to plot or set false

recompute_existing_cluster_stats = true;
stop_early = true;

%% Setup Paths
if hellbender
    addpath 'Raw Data Storage' %#ok<UNRCH>
    addpath 'subfunctions'
    addpath 'clustering subfunctions'
    addpath 'plotting functions'
    output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc');
    patient_IDs = [201907 201908, 201903, 201905, 201906, 201901, 201910, 201915];
else
    addpath('../Z_Raw Data Storage')
    addpath('../subfunctions')
    addpath('../clustering subfunctions')
    addpath('../plotting functions')
    output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
    patient_IDs = [201907]; %#ok<NBRAK>
end

%% Initialize Parallel Pool
if isempty(gcp('nocreate'))
    if hellbender
        parpool('local', 64); %#ok<UNRCH>
    else
        % parpool('local', 2);
    end
end

%% Main Analysis Loop
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    PS_file = get_PS_file(output_folder, patient_ID);
    test_file = get_ES_file(PS_file);
    label_table = PS_file.label_table;
    label_table.anatomical_label = string(label_table.anatomical_label);
    rows_without_nan = get_rows_without_nan(label_table);
    label_table = label_table(rows_without_nan, :);

    all_matrix = squeeze(test_file.all3_ES_matrix(:, :, :, :)); % chan - trials x encoding times x whole trial times  x encoding ID (1-3)
    all_matrix = all_matrix(rows_without_nan, :, :, :);

    % info_file = get_raw_info_file(patient_ID, hellbender);

    for target_enc_ids = 1%:3 % Loop through encoding IDs
        target_matrix = squeeze(all_matrix(:, :, :, target_enc_ids));

        for target_image_ids = 1:2%1:9 % Loop through image IDs
            ctrl_file = get_ctrl_file(PS_file, patient_ID, target_enc_ids, target_image_ids);
            unique_channel_IDs = ctrl_file.unique_channel_IDs;
            num_channels = length(unique_channel_IDs);

            test_table = label_table(ctrl_file.test_rows, :);
            test_matrix = target_matrix(ctrl_file.test_rows, :, :);

            % Pre-allocate results for parallel processing
            channel_results = cell(1, num_channels);

            fprintf("\n computing %d channels for patient%d enc%d image%d\n", num_channels, patient_ID, target_enc_ids, target_image_ids)
            % Pre-slice test_matrix and test_table by channel ID
            pre_sliced_test_matrices = cell(1, num_channels);
            pre_sliced_test_tables = cell(1, num_channels);

            % image_name = info_file.C(target_image_ids);

            for chan_idx = 1:num_channels
                chan_id = unique_channel_IDs(chan_idx);
                chan_test_rows = test_table.channel_ID == chan_id; % Identify rows for this channel
                pre_sliced_test_tables{chan_idx} = test_table(chan_test_rows, :);
                pre_sliced_test_matrices{chan_idx} = test_matrix(chan_test_rows, :, :);
            end
            
            clusters_save_file = strrep(ctrl_file.Properties.Source, 'all_channels_data.mat', 'all_channel_cluster_results.mat');

            if recompute_existing_cluster_stats | ~exist(clusters_save_file, 'file')
                % Parallel loop over channels
                for chan_idx = 2:9:num_channels
                    fprintf('%d ', chan_idx)
                    %clusters_by_anat_local = initialize_clusters_struct(1); % Temporary storage for each channel
                
                    % Access pre-sliced data
                    chan_test_table = pre_sliced_test_tables{chan_idx};
                    chan_test_matrix = pre_sliced_test_matrices{chan_idx};
                
                    % Perform computations for the current channel
                    chan_id = unique_channel_IDs(chan_idx);
                    channel_file = matfile(strrep(ctrl_file.Properties.Source, 'all_channels_data.mat', sprintf('BT_%d', chan_id)));
                    control_matrices = channel_file.control_EMS_matrices(:, :, :, :);
                    control_matrices_std = squeeze(channel_file.control_EMS_matrices_std(:, :, :, :));

                                        %% if num_sig_test_clusters > 0
                    if stop_early
                        plot_test_and_control_EMS(chan_test_matrix, control_matrices, control_matrices_std, chan_id, chan_test_table, target_enc_ids, target_image_ids)
                        % plot_cluster_size_pdf(real_test_cluster_sizes, null_cluster_sizes, nPermutations); % test & null cluster size PDFs
                        % xline(threshold_from_null, 'r', 'LineWidth', 2, 'DisplayName', sprintf('p < %d Threshold from Null', alpha)); % Significance testing
                        % visualize_significant_clusters(real_p_val_matrix, sig_test_cluster_sizes, threshold_from_null, alpha); % view test sig/~sig clusters
                        continue
                    end
                
                    % Process data
                    control_matrices = squeeze(mean(control_matrices, 3));
                    chan_test_matrix = squeeze(mean(chan_test_matrix, 1));
                
                    if size(control_matrices, 2) ~= size(chan_test_matrix, 2)
                        chan_test_matrix = chan_test_matrix(:, 251:640);
                    end
                
                    % Compute clusters
                    [num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, ...
                     null_cluster_sizes, threshold_from_null, real_p_val_matrix] = ...
                     get_clusters(control_matrices, control_matrices_std, chan_test_matrix, nPermutations, alpha, use_z, z_thresh);
    
                    % if num_sig_test_clusters > 0
                        plot_test_and_control_EMS(chan_test_matrix, control_matrices, control_matrices_std, chan_id, chan_test_table, target_enc_ids, target_image_ids)
                        plot_cluster_size_pdf(real_test_cluster_sizes, null_cluster_sizes, nPermutations); % test & null cluster size PDFs
                        xline(threshold_from_null, 'r', 'LineWidth', 2, 'DisplayName', sprintf('p < %d Threshold from Null', alpha)); % Significance testing
                        visualize_significant_clusters(real_p_val_matrix, sig_test_cluster_sizes, threshold_from_null, alpha); % view test sig/~sig clusters
                    % end
                
                    % Store Results
                    anat = unique(chan_test_table.anatomical_label);
                    if numel(anat) ~= 1
                        error("Each test data should have a single anatomical label.");
                    end
                
                    clusters_by_anat_local = store_cluster_results(clusters_by_anat_local, 1, anat, ...
                        num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes);
                
                    % Save results to the cell array
                    channel_results{chan_idx} = clusters_by_anat_local;
                end
                save(clusters_save_file, "channel_results", "-v7.3")
            else
                data = load(clusters_save_file); %#ok<UNRCH>
                channel_results = data.channel_results;
                clear data
            end

            % clusters = fix_structure(channel_results);
            % clusters.anat = string(clusters.anat);
            % clusters.num_sig_test_clusters = cell2mat(clusters.num_sig_test_clusters);
            % clusters.sig_test_cluster_sizes = cell2mat(clusters.sig_test_cluster_sizes);
            % clusters.real_test_cluster_sizes = cell2mat(clusters.real_test_cluster_sizes);

            % [clusters_by_uniq_anat] = sort_by_uniq_anat(clusters);
            

            % plot_cluster_sizes_by_uniq_anat(clusters_by_uniq_anat)

        end
    end
end

%% Function Definitions
function clusters = initialize_clusters_struct(num_channels)
    clusters.anat = cell(num_channels, 1);
    clusters.num_sig_test_clusters = cell(num_channels, 1);
    clusters.sig_test_cluster_sizes = cell(num_channels, 1);
    clusters.real_test_cluster_sizes = cell(num_channels, 1);
end

function clusters = store_cluster_results(clusters, idx, anat, num_clusters, sig_sizes, real_sizes)
    clusters.anat{idx} = anat;
    clusters.num_sig_test_clusters{idx} = num_clusters;
    clusters.sig_test_cluster_sizes{idx} = sig_sizes;
    clusters.real_test_cluster_sizes{idx} = real_sizes;
end

function clusters_main = merge_cluster_results(clusters_main, clusters_local)
    for idx = 1:length(clusters_local.anat)
        if ~isempty(clusters_local.anat{idx})
            clusters_main.anat{end + 1} = clusters_local.anat{idx};
            clusters_main.num_sig_test_clusters{end + 1} = clusters_local.num_sig_test_clusters{idx};
            clusters_main.sig_test_cluster_sizes{end + 1} = clusters_local.sig_test_cluster_sizes{idx};
            clusters_main.real_test_cluster_sizes{end + 1} = clusters_local.real_test_cluster_sizes{idx};
        end
    end
end

function [num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, null_cluster_sizes, threshold_from_null, real_p_val_matrix] = get_clusters(control_matrices, control_matrices_std, test_matrix, nPermutations, alpha, use_z, z_thresh)
    sanity_check_matrices(test_matrix, control_matrices);

    % Compute real test's z-scored or p-value matrices
    if use_z
        real_p_val_matrix = [];
        % test_z_matrix = get_z_mat(test_matrix, control_matrices);
        test_z_matrix = z_test_by_control(test_matrix, control_matrices, control_matrices_std);
        real_test_cluster_sizes = compute_cluster_sizes_z(test_z_matrix, z_thresh);
    else
        % test_z_matrix = [];
        real_p_val_matrix = get_p_mat(test_matrix, control_matrices);
        real_test_cluster_sizes = compute_cluster_sizes(real_p_val_matrix, alpha);
    end

    % Concatenate control and test matrices for null distribution permutations
    all_avg_matrices = cat(3, control_matrices, test_matrix);

    % Run permutations
    null_cluster_sizes = run_permutations(nPermutations, all_avg_matrices, alpha, use_z, z_thresh);

    % Concatenate null cluster sizes
    null_cluster_sizes = concatenate_cluster_sizes(null_cluster_sizes);

    % Find significant clusters
    threshold_from_null = prctile(null_cluster_sizes, 100 * (1 - alpha));
    sig_test_cluster_sizes = real_test_cluster_sizes(real_test_cluster_sizes > threshold_from_null);
    num_sig_test_clusters = numel(sig_test_cluster_sizes);
end

function p_values = get_p_mat(test_mat, control_mats)
    % Pre-allocate p_values
    % p_values = zeros(size(test_mat));
    
    % Reshape test_mat to broadcast over the third dimension
    test_mat_expanded = reshape(test_mat, size(test_mat, 1), size(test_mat, 2), 1);
    
    % Vectorized percentile computation
    p_values = sum(control_mats > test_mat_expanded, 3) ./ size(control_mats, 3);
end

function z_values = get_z_mat(test_mat, control_mats)
    % Pre-allocate z_values
    % z_values = zeros(size(test_mat));
    
    % Compute mean and std along the third dimension in one step
    dist_mean = mean(control_mats, 3);
    dist_std = std(control_mats, 0, 3); % 0 specifies normalization by N-1
    
    % Vectorized z-score computation
    z_values = (test_mat - dist_mean) ./ dist_std;
end

%alternative if provided the original std matrix
function test_z = z_test_by_control(test_EMS_avg, control_EMS_avg, control_EMS_std)
    control_EMS_avg = squeeze(mean(control_EMS_avg, 3)); % average EMS avg across samples
    control_EMS_std = squeeze(mean(control_EMS_std, 3)); % average EMS std avg across samples
    test_z = (test_EMS_avg - control_EMS_avg) ./ control_EMS_std;
end

function clusterSizesNull = run_permutations(nPermutations, all_avg_matrices, alpha, use_z, z_thresh)
    clusterSizesNull = cell(nPermutations, 1);  % Preallocate

    for i = 1:nPermutations
        [test_mat, control_mats] = pick_rand_test(all_avg_matrices);

        if use_z
            % Compute z-scored matrix
            z_scored_mat = get_z_mat(test_mat, control_mats);
            clusterSizes = compute_cluster_sizes_z(z_scored_mat, z_thresh);
        else
            % Compute p-value matrix
            p_val_matrix = get_p_mat(test_mat, control_mats);
            clusterSizes = compute_cluster_sizes(p_val_matrix, alpha);
        end

        clusterSizesNull{i} = clusterSizes;
    end
end

function clusters = fix_structure(channel_results)
    % converts cell array of structs to struct of cell arrays
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
end
