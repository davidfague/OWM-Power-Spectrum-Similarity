%% Initialization Parameters
params = struct();
params.hellbender = false;
params.within_item_btwn_trial = true; % false: within_item_within_trial test vs between_item_between_trial; true: within_item_btwn_trial test vs btwn_item_btwn_trial
params.nPermutations = 1000;
params.use_z = true; % false uses p
params.alpha = 0.05;
params.z_thresh = 1.96;
params.recompute_existing_cluster_stats = false;
params.debug = false; % Set to false to suppress debug statements
params.output_folder = get_output_folder(params.hellbender);
params.patient_IDs = get_patient_IDs(params.hellbender);

%% Setup Paths and Parallel Pool
setup_paths(params.hellbender);
initialize_parallel_pool(params.hellbender);

%% Main Analysis Loop
results_by_anat = analyze_patients(params);

%% Plot Results
plot_results_by_anatomy(results_by_anat);

%% Function Definitions
function output_folder = get_output_folder(hellbender)
    if hellbender
        output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc');
    else
        output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
    end
end

function patient_IDs = get_patient_IDs(hellbender)
    if hellbender
        patient_IDs = [201907, 201908, 201903, 201905, 201906, 201901, 201910, 201915];
    else
        patient_IDs = [201907];
    end
end

function setup_paths(hellbender)
    if hellbender
        addpath 'Raw Data Storage';
        addpath 'subfunctions';
        addpath 'clustering subfunctions';
    else
        addpath('../Z_Raw Data Storage');
        addpath('../subfunctions');
        addpath('../clustering subfunctions');
    end
end

function initialize_parallel_pool(hellbender)
    if isempty(gcp('nocreate'))
        if hellbender
            parpool('local', 64);
        else
            parpool('local', 1);
        end
    end
end

function clusters_by_anat = analyze_patients(params)
    % Initialize structure for results
    clusters_by_anat = struct();

    for patient_ID = params.patient_IDs
        fprintf('Processing Patient %d...\n', patient_ID);
        [PS_file, label_table] = load_patient_files(params.output_folder, patient_ID);

        for enc_ID = 1:3
            fprintf('  Encoding ID: %d\n', enc_ID);
            for image_ID = 1:9
                fprintf('    Image ID: %d\n', image_ID);

                % Dynamically load the test file
                test_file = get_test_file(PS_file, enc_ID, image_ID, params.within_item_btwn_trial);

                % Compute or load clusters
                clusters_save_file = get_clusters_save_file(PS_file, enc_ID, image_ID);
                recompute_needed = params.recompute_existing_cluster_stats || ~exist(clusters_save_file, 'file');
                results_missing = false;

                if ~recompute_needed
                    try
                        data = load(clusters_save_file, 'results');
                        if isfield(data, 'results')
                            patient_clusters = data.results;
                        else
                            warning("Saved file '%s' does not contain 'results'. Recomputing...", clusters_save_file);
                            results_missing = true;
                        end
                    catch
                        warning("Failed to load '%s'. Recomputing...", clusters_save_file);
                        recompute_needed = true;
                    end
                end

                if recompute_needed || results_missing
                    fprintf("Recomputing clusters for enc_ID: %d, image_ID: %d\n", enc_ID, image_ID);
                    patient_clusters = compute_clusters_for_image(params, PS_file, test_file, label_table, enc_ID, image_ID);
                    save(clusters_save_file, 'results', '-v7.3');
                else
                    fprintf("Loading existing clusters for enc_ID: %d, image_ID: %d\n", enc_ID, image_ID);
                end

                % Merge results into the main structure
                clusters_by_anat = merge_clusters(clusters_by_anat, patient_clusters);
            end
        end
    end
end

function [PS_file, label_table] = load_patient_files(output_folder, patient_ID)
    PS_file = get_PS_file(output_folder, patient_ID);
    label_table = PS_file.label_table;
    label_table.anatomical_label = string(label_table.anatomical_label);
    label_table = label_table(get_rows_without_nan(label_table), :);
end

% function results = process_image(params, PS_file, test_file, label_table, enc_ID, image_ID)
%     % Load control matrix (btwn_item_btwn_trial remains the same)
%     ctrl_file = get_ctrl_file(PS_file, params.patient_IDs, enc_ID, image_ID);
%     preloaded_ctrl_matrices = preload_control_matrices(ctrl_file, ctrl_file.unique_channel_IDs);
% 
%     % Pre-slice test data
%     unique_channel_IDs = ctrl_file.unique_channel_IDs;
%     if params.within_item_btwn_trial
%         pre_sliced_data = pre_slice_test_data(label_table, test_file.control_EMS_matrices, unique_channel_IDs, params.debug);
%     else
%         pre_sliced_data = pre_slice_test_data(label_table, squeeze(test_file.all3_ES_matrix(:, :, :, enc_ID)), unique_channel_IDs, params.debug);
%     end
% 
%     % Perform cluster analysis
%     clusters_save_file = strrep(ctrl_file.Properties.Source, 'all_channels_data.mat', 'all_channel_cluster_results.mat');
%     recompute_needed = params.recompute_existing_cluster_stats || ~exist(clusters_save_file, 'file');
%     results_missing = false;
% 
%     if ~recompute_needed
%         try
%             data = load(clusters_save_file, 'results');
%             if isfield(data, 'results')
%                 results = data.results;
%             else
%                 warning("Saved file '%s' does not contain 'results'. Recomputing...", clusters_save_file);
%                 results_missing = true;
%             end
%         catch
%             warning("Failed to load '%s'. Recomputing...", clusters_save_file);
%             recompute_needed = true;
%         end
%     end
% 
%     if recompute_needed || results_missing
%         fprintf("Recomputing clusters for enc_ID: %d, image_ID: %d\n", enc_ID, image_ID);
%         results = compute_clusters(pre_sliced_data, preloaded_ctrl_matrices, params);
%         save(clusters_save_file, "results", "-v7.3");
%     else
%         fprintf("Loading existing clusters for enc_ID: %d, image_ID: %d\n", enc_ID, image_ID);
%     end
% end

function pre_sliced_data = pre_slice_test_data(label_table, test_file, unique_channel_IDs, enc_ID, params)
    % Pre-slice test data by channel
    num_channels = numel(unique_channel_IDs);

    % Initialize storage
    pre_sliced_data = struct();
    pre_sliced_data.tables = cell(1, num_channels);
    pre_sliced_data.matrices = cell(1, num_channels);

    for chan_idx = 1:num_channels
        chan_id = unique_channel_IDs(chan_idx);
        chan_rows = label_table.channel_ID == chan_id;

        % Handle within-item/between-item test file structure
        if params.within_item_btwn_trial
            % For within-item trials, load data per channel
            channel_file_path = strrep(test_file.Properties.Source, ...
                                       'all_channels_data.mat', sprintf('BT_%d.mat', chan_id));
            if ~exist(channel_file_path, 'file')
                warning("Channel-specific file '%s' does not exist. Skipping channel ID %d.", channel_file_path, chan_id);
                pre_sliced_data.matrices{chan_idx} = [];
            else
                channel_file = matfile(channel_file_path);
                if any(strcmp('control_EMS_matrices', who(channel_file)))
                    pre_sliced_data.matrices{chan_idx} = squeeze(channel_file.control_EMS_matrices(:, :, :, :));
                else
                    error("'control_EMS_matrices' not found in '%s'.", channel_file_path);
                end
            end
        else
            % For standard encoding/maintenance trials, use test_file directly
            if isfield(test_file, 'all3_ES_matrix')
                pre_sliced_data.matrices{chan_idx} = squeeze(test_file.all3_ES_matrix(chan_rows, :, :, enc_ID));
            else
                error("'all3_ES_matrix' not found in the provided test file.");
            end
        end

        % Subset label_table rows
        pre_sliced_data.tables{chan_idx} = label_table(chan_rows, :);

        if params.debug
            fprintf("Channel %d: Rows=%d\n", chan_idx, nnz(chan_rows));
        end
    end
end

% function aggregated_results = compute_clusters(pre_sliced_data, preloaded_ctrl_matrices, params)
%     num_channels = length(pre_sliced_data.tables);
% 
%     % Initialize result structure
%     aggregated_results = struct();
%     aggregated_results.anat = [];
%     aggregated_results.num_sig_test_clusters = [];
%     aggregated_results.sig_test_cluster_sizes = [];
%     aggregated_results.real_test_cluster_sizes = [];
% 
%     parfor chan_idx = 1:num_channels
%         % Perform cluster computation
%         chan_data = pre_sliced_data.matrices{chan_idx};
%         control_matrices = preloaded_ctrl_matrices{chan_idx};
% 
%         if isempty(control_matrices)
%             warning("Skipping channel index %d due to missing control matrices.", chan_idx);
%             continue;
%         end
% 
%         [clusters, anat] = compute_channel_clusters(chan_data, control_matrices, params);
% 
%         % Combine results into a temporary struct
%         temp_results = struct();
%         temp_results.anat = anat;
%         temp_results.num_sig_test_clusters = clusters.num_sig_test_clusters;
%         temp_results.sig_test_cluster_sizes = clusters.sig_test_cluster_sizes;
%         temp_results.real_test_cluster_sizes = clusters.real_test_cluster_sizes;
% 
%         % Update aggregated results (using parallel-safe accumulation)
%         aggregated_results = merge_results(aggregated_results, temp_results);
%     end
% end

% function [clusters, anat] = compute_channel_clusters(test_matrix, control_matrix, params)
%     % Compute clusters
%     clusters = compute_clusters_from_data(control_matrix, test_matrix, params);
% 
%     % Extract anatomy information
%     anat = unique(params.label_table.anatomical_label); % Assuming label_table is passed in params
%     assert(isscalar(anat), "Each test data must have a single anatomical label.");
% end

% function clusters = compute_clusters_from_data(pre_sliced_test_tables, pre_sliced_test_matrices, unique_channel_IDs, ctrl_file, params)
%     % Computes clusters for test and control matrices, storing results directly into a structured format.
%     %
%     % Inputs:
%     % - pre_sliced_test_tables: Cell array of channel-specific test data tables
%     % - pre_sliced_test_matrices: Cell array of channel-specific test data matrices
%     % - unique_channel_IDs: List of unique channel IDs
%     % - ctrl_file: Control file with precomputed control matrices
%     % - params: Parameters for clustering and statistical analysis
%     %
%     % Output:
%     % - clusters: Struct containing results organized by anatomy, patient, encoding, and image
% 
%     % Initialize clusters structure
%     clusters = struct();
% 
%     % File for saving results
%     clusters_save_file = strrep(ctrl_file.Properties.Source, 'all_channels_data.mat', 'all_channel_cluster_results.mat');
% 
%     if params.recompute_existing_cluster_stats || ~exist(clusters_save_file, 'file')
%         % Parallel loop over channels
%         parfor chan_idx = 1:numel(unique_channel_IDs)
%             fprintf('%d ', chan_idx);
% 
%             % Access pre-sliced test data for this channel
%             chan_test_table = pre_sliced_test_tables{chan_idx};
%             chan_test_matrix = pre_sliced_test_matrices{chan_idx};
% 
%             % Get control matrices for this channel
%             chan_id = unique_channel_IDs(chan_idx);
%             channel_file = matfile(strrep(ctrl_file.Properties.Source, 'all_channels_data.mat', sprintf('BT_%d', chan_id)));
%             control_matrices = channel_file.control_EMS_matrices(:, :, :, :);
%             control_matrices_std = channel_file.control_EMS_matrices_std(:, :, :, :);
% 
%             % Process data
%             control_matrices = squeeze(mean(control_matrices, 3));
%             chan_test_matrix = squeeze(mean(chan_test_matrix, 1));
% 
%             % Ensure matching dimensions
%             if size(control_matrices, 2) ~= size(chan_test_matrix, 2)
%                 chan_test_matrix = chan_test_matrix(:, 251:640);
%             end
% 
%             % Compute clusters
%             [num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, null_cluster_sizes, threshold_from_null] = ...
%                 get_clusters(control_matrices, control_matrices_std, chan_test_matrix, params);
% 
%             % Organize results by anatomy, patient, encoding, and image
%             anat = unique(chan_test_table.anatomical_label);
%             if numel(anat) ~= 1
%                 error("Each test data should have a single anatomical label.");
%             end
% 
%             patient = sprintf('Patient_%d', params.patient_IDs(chan_idx));
%             encoding = sprintf('Encoding_%d', params.enc_ID);
%             image = sprintf('Image_%d', params.image_ID);
% 
%             % Store results in the main structure
%             clusters = store_cluster_results(clusters, patient, anat, encoding, image, ...
%                 num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, null_cluster_sizes, threshold_from_null);
%         end
% 
%         % Save results
%         save(clusters_save_file, "clusters", "-v7.3");
%     else
%         % Load precomputed results
%         data = load(clusters_save_file);
%         clusters = data.clusters;
%         clear data;
%     end
% end

% function aggregated_results = merge_results(aggregated_results, temp_results)
%     % Merge individual channel results into the aggregated struct
%     aggregated_results.anat = [aggregated_results.anat; temp_results.anat];
%     aggregated_results.num_sig_test_clusters = [aggregated_results.num_sig_test_clusters; temp_results.num_sig_test_clusters];
%     aggregated_results.sig_test_cluster_sizes = [aggregated_results.sig_test_cluster_sizes; temp_results.sig_test_cluster_sizes];
%     aggregated_results.real_test_cluster_sizes = [aggregated_results.real_test_cluster_sizes; temp_results.real_test_cluster_sizes];
% end

function plot_results_by_anatomy(results_by_anat)
    unique_anatomies = fieldnames(results_by_anat);
    all_sizes = [];
    groups = [];

    for i = 1:numel(unique_anatomies)
        anat = unique_anatomies{i};
        cluster_sizes = results_by_anat.(anat).real_test_cluster_sizes;
        all_sizes = [all_sizes; cluster_sizes];
        groups = [groups; repmat({anat}, length(cluster_sizes), 1)];
    end

    % Generate boxplot
    figure;
    boxplot(all_sizes, groups, 'LabelOrientation', 'inline');
    xlabel('Anatomical Labels');
    ylabel('Cluster Size');
    title('Cluster Sizes Grouped by Anatomical Regions');
    set(gca, 'XTickLabelRotation', 45);
    grid on;
end

function [preloaded_ctrl_matrices, preloaded_ctrl_matrices_std] = preload_control_matrices(ctrl_file, unique_channel_IDs)
    % Pre-slice control matrices by channel
    num_channels = length(unique_channel_IDs);
    preloaded_ctrl_matrices = cell(1, num_channels);
    preloaded_ctrl_matrices_std = cell(1, num_channels);

    for chan_idx = 1:num_channels
        chan_id = unique_channel_IDs(chan_idx);

        % Construct the channel-specific file path
        channel_file_path = strrep(ctrl_file.Properties.Source, ...
                                   'all_channels_data.mat', sprintf('BT_%d.mat', chan_id));
        if ~exist(channel_file_path, 'file')
            warning("Channel file '%s' does not exist. Skipping channel ID %d.", channel_file_path, chan_id);
            preloaded_ctrl_matrices{chan_idx} = [];
            preloaded_ctrl_matrices_std{chan_idx} = [];
            continue;
        end

        % Load the control matrices and their standard deviations
        try
            channel_file = matfile(channel_file_path);
            preloaded_ctrl_matrices{chan_idx} = channel_file.control_EMS_matrices(:, :, :, :);
            preloaded_ctrl_matrices_std{chan_idx} = channel_file.control_EMS_matrices_std(:, :, :, :);
        catch ME
            warning("Failed to load control matrices for channel ID %d: %s", chan_id, ME.message);
            preloaded_ctrl_matrices{chan_idx} = [];
            preloaded_ctrl_matrices_std{chan_idx} = [];
        end
    end
end


function test_file = get_test_file(PS_file, enc_ID, image_ID, within_item_vs_btwn)
    % Determine test file path based on within_item_vs_btwn
    if within_item_vs_btwn %true: within_item_btwn_trial test vs btwn_item_btwn_trial
        test_file_path = fullfile(strrep(PS_file.Properties.Source, '.mat', ...
            sprintf('/corrWitemMaint vs corrWitemEnc/between_trials_enc%d_image%d/all_channels_data.mat', enc_ID, image_ID)));

        % Check if the file exists
        if ~exist(test_file_path, 'file')
            
            error("Test file not found: %s", test_file_path);
        end
    
        % Load the test file
        test_file = matfile(test_file_path);

    else % within_item_within_trial test vs between_item_between_trial
        test_file = get_ES_file(PS_file); % returns a matfile
    end
end

function clusters_by_anat = compute_clusters_for_image(params, PS_file, test_file, label_table, enc_ID, image_ID)
    % Prepare control file and pre-sliced data
    ctrl_file = get_ctrl_file(PS_file, params.patient_IDs, enc_ID, image_ID);
    [preloaded_ctrl_matrices, preloaded_ctrl_matrices_std] = preload_control_matrices(ctrl_file, ctrl_file.unique_channel_IDs);

    % Pre-slice test data
    unique_channel_IDs = ctrl_file.unique_channel_IDs;
    pre_sliced_data = pre_slice_test_data(label_table, test_file, unique_channel_IDs, enc_ID, params);

    % Temporary storage for results
    num_channels = numel(unique_channel_IDs);
    clusters_temp = cell(1, num_channels);

    % Parallel loop over channels
    for chan_idx = 1:num_channels
        fprintf('%d ', chan_idx);

        % Extract data for the current channel
        chan_test_data = pre_sliced_data.matrices{chan_idx};
        chan_control_matrices = preloaded_ctrl_matrices{chan_idx};
        chan_control_matrices_std = preloaded_ctrl_matrices_std{chan_idx};

        if isempty(chan_control_matrices) || isempty(chan_control_matrices_std)
            warning("Skipping channel index %d due to missing control matrices or std.", chan_idx);
            continue;
        end

        % Compute clusters for the channel
        [num_sig_clusters, sig_cluster_sizes, real_cluster_sizes, null_cluster_sizes, threshold] = ...
            get_clusters(chan_control_matrices, chan_control_matrices_std, chan_test_data, params);

        % Extract anatomy information
        anat = unique(pre_sliced_data.tables{chan_idx}.anatomical_label);
        if numel(anat) ~= 1
            error("Each test data must have a single anatomical label.");
        end

        % Store results for this channel
        clusters_temp{chan_idx} = struct(...
            'anat', anat, ...
            'num_sig_clusters', num_sig_clusters, ...
            'sig_cluster_sizes', sig_cluster_sizes, ...
            'real_cluster_sizes', real_cluster_sizes, ...
            'null_cluster_sizes', null_cluster_sizes, ...
            'threshold', threshold, ...
            'enc_ID', enc_ID, ...
            'image_ID', image_ID ...
        );
    end

    % Aggregate results after the `parfor` loop
    clusters_by_anat = aggregate_cluster_results(clusters_temp);
end

% function clusters = initialize_clusters_struct()
%     % Initialize an empty structure for clusters
%     clusters = struct();
% end

% function clusters = store_cluster_results(clusters, anat, num_sig_clusters, sig_cluster_sizes, real_cluster_sizes, null_cluster_sizes, threshold, enc_ID, image_ID)
%     % Ensure the anatomy field exists
%     if ~isfield(clusters, anat)
%         clusters.(anat) = struct();
%     end
% 
%     % Store results under the appropriate encoding and image
%     field_name = sprintf('Enc%d_Image%d', enc_ID, image_ID);
%     clusters.(anat).(field_name) = struct(...
%         'num_sig_clusters', num_sig_clusters, ...
%         'sig_cluster_sizes', sig_cluster_sizes, ...
%         'real_cluster_sizes', real_cluster_sizes, ...
%         'null_cluster_sizes', null_cluster_sizes, ...
%         'threshold', threshold ...
%     );
% end

function clusters_main = merge_clusters(clusters_main, clusters_local)
    % Merge local clusters into the main structure
    anat_fields = fieldnames(clusters_local);
    for i = 1:numel(anat_fields)
        anat = anat_fields{i};
        if ~isfield(clusters_main, anat)
            clusters_main.(anat) = clusters_local.(anat);
        else
            main_fields = fieldnames(clusters_local.(anat));
            for j = 1:numel(main_fields)
                field_name = main_fields{j};
                clusters_main.(anat).(field_name) = clusters_local.(anat).(field_name);
            end
        end
    end
end

function [num_sig_clusters, sig_cluster_sizes, real_cluster_sizes, null_cluster_sizes, threshold] = get_clusters(control_matrices, control_matrices_std, test_matrix, params)
    % Compute clusters for the provided matrices
    if params.use_z
        test_z_matrix = z_test_by_control(test_matrix, control_matrices, control_matrices_std);
        real_cluster_sizes = compute_cluster_sizes_z(test_z_matrix, params.z_thresh);
    else
        real_p_val_matrix = get_p_mat(test_matrix, control_matrices);
        real_cluster_sizes = compute_cluster_sizes(real_p_val_matrix, params.alpha);
    end

    % Permutation test for null distribution
    null_cluster_sizes = run_permutations(params, cat(3, squeeze(control_matrices), test_matrix));

    % Significance threshold and significant clusters
    flat_cluster_sizes = get_flat_cluster_sizes(null_cluster_sizes);
    threshold = prctile(flat_cluster_sizes, 100 * (1 - params.alpha));
    sig_cluster_sizes = real_cluster_sizes(real_cluster_sizes > threshold);
    num_sig_clusters = numel(sig_cluster_sizes);
end

function flat_cluster_sizes = get_flat_cluster_sizes(null_cluster_sizes)
    % Reshape each non-empty cell element into a row vector
    reshaped_sizes = cellfun(@(x) reshape(x, 1, []), null_cluster_sizes, 'UniformOutput', false);
    
    % Remove empty arrays
    reshaped_sizes = reshaped_sizes(~cellfun('isempty', reshaped_sizes));
    
    % Concatenate all reshaped arrays into one numerical array
    flat_cluster_sizes = [reshaped_sizes{:}];
    
    % Check if the resulting array is empty
    if isempty(flat_cluster_sizes)
        error('No valid cluster sizes found in null_cluster_sizes.');
    end
end

function clusters = aggregate_cluster_results(clusters_temp)
    % Initialize an empty structure for aggregated results
    clusters = struct();

    for idx = 1:numel(clusters_temp)
        if isempty(clusters_temp{idx})
            continue;
        end

        % Extract data for this channel
        chan_data = clusters_temp{idx};
        anat = chan_data.anat;
        field_name = sprintf('Enc%d_Image%d', chan_data.enc_ID, chan_data.image_ID);

        % Ensure anatomy field exists
        if ~isfield(clusters, anat)
            clusters.(anat) = struct();
        end

        % Store results under the appropriate encoding and image
        clusters.(anat).(field_name) = rmfield(chan_data, {'anat', 'enc_ID', 'image_ID'});
    end
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

function clusterSizesNull = run_permutations(params, all_avg_matrices)
    clusterSizesNull = cell(params.nPermutations, 1);  % Preallocate

    for i = 1:params.nPermutations
        [test_mat, control_mats] = pick_rand_test(all_avg_matrices);

        if params.use_z
            % Compute z-scored matrix
            z_scored_mat = get_z_mat(test_mat, control_mats);
            clusterSizes = compute_cluster_sizes_z(z_scored_mat, params.z_thresh);
        else
            % Compute p-value matrix
            p_val_matrix = get_p_mat(test_mat, control_mats);
            clusterSizes = compute_cluster_sizes(p_val_matrix, params.alpha);
        end

        clusterSizesNull{i} = clusterSizes;
    end
end

function clusters_save_file = get_clusters_save_file(PS_file, enc_ID, image_ID)
    clusters_save_file = strrep(PS_file.Properties.Source, ...
        'all_channels_data.mat', ...
        sprintf('Enc%d_Image%d_cluster_results.mat', enc_ID, image_ID));
end
