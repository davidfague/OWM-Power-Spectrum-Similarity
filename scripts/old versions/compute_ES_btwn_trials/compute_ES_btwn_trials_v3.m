%% compute control: correct-item-encoding matched with random non-item, all3 correct maintenance
%% v3: use get_parameters.m
%% v2: remove the without_nan subsetting of label_table and PSVs to keep indexing inconsistent; also begin saving only unique results
% also only save the unique trial pairs, do not oversample.

params = get_parameters();
    % temp params
params.within_item = true; % computes control as either within item or between items

if params.hellbender
    % addpath 'Raw Data Storage'               %#ok<UNRCH>
    addpath 'subfunctions'
else
    % addpath('../Z_Raw Data Storage')      %#ok<UNRCH>    
    addpath('../subfunctions')
end

%% main loop
for idx = 1:length(params.patient_IDs)
    patient_ID = params.patient_IDs(idx);
    PS_file = get_PS_file(params.output_folder, patient_ID);

    % % % compute correlations
    for target_enc_ids = 1%:3 % can be 1:3
        for target_image_ids = [1,3] % can be 1:9
            fprintf("computing patient%d, enc%d, image%d \n", patient_ID, target_enc_ids, target_image_ids)
            [all_channels_save_file] = compute_all_between_trial_similarities(PS_file, target_enc_ids, target_image_ids, params);
        end
    end

end

%% main computations
function [all_channels_save_file] = compute_all_between_trial_similarities(PS_file, target_enc_ids, target_image_ids, params)
% would need reworking to accept all 3 encodings at once. Have to select
% enc_win_IDs based on the encID where the thing is found.

    % Select encoding window IDs based on target encoding IDs
    switch target_enc_ids
        case 1
            enc_win_IDs = params.enc1_win_IDs;
        case 2
            enc_win_IDs = params.enc2_win_IDs;
        case 3
            enc_win_IDs = params.enc3_win_IDs;
    end

    % data_file.all_windowed_mean_PS_vectors: size = all_window_IDs, frequencies, channels x trials 
    % label_table = data_file.label_table;

    % Get the cleaned data from the data file
    [label_table, all_window_mean_PS_vectors] = get_cleaned_data(PS_file);
    %% testing
    % target_trial = label_table.trial_ID == 43;
    % target_channel = label_table.channel_ID==3;
    % target_row = target_channel & target_trial;
    % sum(isnan(all_window_mean_PS_vectors(:,:,target_row)))
    %%

    % Pool the target image and encoding correctness for selected encIDs
    test_rows = filter_rows(label_table, target_image_ids, target_enc_ids);
    item_cor_table = label_table(test_rows, :);
    item_cor_PSVs = all_window_mean_PS_vectors(:, :, test_rows);

    % Get control data (all 3 encoding correct without the target image)
    if params.within_item
        control_rows = test_rows;
    else
        control_rows = filter_all_correct_without_target_image(label_table, target_image_ids);
    end
    control_table = label_table(control_rows, :);
    control_PSVs = all_window_mean_PS_vectors(:, :, control_rows);
    clearvars all_window_mean_PS_vectors label_table

    % Identify unique channels
    unique_channel_IDs = unique(item_cor_table.channel_ID);

    % Set up save folder
    if params.within_item
        save_folder = fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('/corrWitemMaint vs corrWitemEnc/between_trials_enc%d_image%d', target_enc_ids, target_image_ids)));
    else
        save_folder = fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('/corrWOitemMaint vs corrcWitemEnc/between_trials_enc%d_image%d', target_enc_ids, target_image_ids)));
    end
    % save_folder = strrep(PS_file.Properties.Source, '.mat', '_all3corrWOitem1_and_enc1correctWitem1');
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    % Initialize storage structures for all channels
    % all_test_EMS_matrices_by_chan = cell(length(unique_channel_IDs), 1); % no longer recomputing test, gathering it in next script
    % all_control_EMS_matrices_by_chan = cell(length(unique_channel_IDs), 1);

    %% Loop through each unique channel
    for chan_idx = 1:length(unique_channel_IDs)
        chan_id = unique_channel_IDs(chan_idx);
        save_file = fullfile(save_folder, sprintf('BT_%d.mat', chan_id));
        
        % Filter data for the current channel ID
        chan_test_rows = filter_by_channel_id(item_cor_table, chan_id);
        % [chan_test_table, chan_test_PSVs] = filter_table_PSVs(item_cor_table, item_cor_PSVs, chan_test_rows);
        [~, chan_test_PSVs] = filter_table_PSVs(item_cor_table, item_cor_PSVs, chan_test_rows);

        chan_ctrl_rows = filter_by_channel_id(control_table, chan_id);
        % [chan_ctrl_table, chan_ctrl_PSVs] = filter_table_PSVs(control_table, control_PSVs, chan_ctrl_rows);
        [~, chan_ctrl_PSVs] = filter_table_PSVs(control_table, control_PSVs, chan_ctrl_rows);
        %% check nans
        % nanPresence = any(isnan(chan_ctrl_PSVs), 2);
        % figure;
        % imagesc(squeeze(nanPresence));
        % colormap([1 1 1; 0 0 0]);
        %%
        clear chan_ctrl_rows

        % % Compute test matrices
        num_test_trials = size(chan_test_PSVs, 3);
        % test_EMS_matrices = zeros(length(enc_win_IDs), length(maint_win_IDs), num_trials);
        % for j = 1:num_trials
        %     % Compute test matrix for each trial
        %     test_EMS_matrices(:, :, j) = compute_similarity_matrix(chan_test_PSVs(:, :, j), enc_win_IDs, maint_win_IDs);
        % end

        % Compute control matrices with random pairs
        numSamples = 1000; % increasing from 100 to 1000
        num_ctrl_trials = size(chan_ctrl_PSVs, 3);
                % num_possible_unique_pairs = num_test_trials * num_ctrl_trials;% 20 * ~75 = ~1500
        % num_pairs = num_test_trials * numSamples % = 20 * 1000 = 20000

        %% testing computing only unique trial pairs then duplicating to get entire matrix.

        rand_trials = randi(num_ctrl_trials, num_test_trials, numSamples); % num_test_trials x numSamples with ctrl_trial values
        % Reshape rand_trials to obtain a list of (test_trial, control_trial) pairs
        pairs_list = [repmat((1:num_test_trials)',numSamples, 1), rand_trials(:)];

        if params.within_item % prevent matching trial with itself for within-item.
            pairs_list = filter_trial_pairs_for_within_item(pairs_list, num_ctrl_trials);
        end

        [unique_pairs, ~, idx_map] = unique(pairs_list, 'rows'); % Identify unique pairs and get indices mapping each pair to its unique entry
        unique_results = nan(length(enc_win_IDs), length(params.maint_win_IDs), size(unique_pairs, 1)); % Initialize a temporary array to store computations for each unique pair
        clear rand_trials pairs_list num_ctrl_trials
        %% memory management
        % sum([whos().bytes]) / (1024^3) = 3.1 GB
        % vars_before = {whos().name};
        %clearvars -except vars_before unique_results chan_test_PSVs chan_ctrl_PSVs enc_win_IDs maint_win_IDs unique_pairs numSamples num_test_trials chan_idx all_control_EMS_matrices_by_chan save_folder target_image_ids target_enc_ids test_rows control_rows unique_channel_IDs idx_map item_cor_table item_cor_PSVs control_table control_PSVs
        % vars_after = {whos().name};
        % deleted_vars = setdiff(vars_before, vars_after);
        % disp('Deleted Variables')
        % disp(deleted_vars)
        %% Compute similarity matrices for each unique pair
        fprintf("computing %d unique pairs instead of all %d*%d for chan%d\n", size(unique_pairs, 1), num_test_trials,numSamples, chan_id)
        for k = 1:size(unique_pairs, 1)
            j = unique_pairs(k, 1);        % test trial index
            randtrial = unique_pairs(k, 2); % control trial index

            if any(any(isnan(chan_test_PSVs(:, :, j)))) || any(any(isnan(chan_ctrl_PSVs(:, :, randtrial)))) || all(all(chan_ctrl_PSVs(:, :, randtrial)==0))|| all(all(chan_test_PSVs(:, :, j)==0))
                % error("nan detected")
                unique_results(:, :, k) = nan(size(unique_results, [1 2]));
            end
            
            % Compute and store the result for the unique pair
            unique_results(:, :, k) = compute_BT_similarity_matrix(chan_test_PSVs(:, :, j), chan_ctrl_PSVs(:, :, randtrial), enc_win_IDs, params.maint_win_IDs);
        end
        clearvars j randtrial chan_ctrl_PSVs chan_test_PSVs
        %% new

        control_EMS_matrices = unique_results;
        control_EMS_matrices_std = std(control_EMS_matrices, 0, 3);
        % save(save_file, "chan_id", "control_EMS_matrices", "control_EMS_matrices_std", "chan_test_rows", "-v7.3");
        save(save_file, "chan_id", "control_EMS_matrices", "control_EMS_matrices_std", "chan_test_rows");
        %% estimating memory usage per # parallel workers
        % Example sizes for variables on each worker
        % memory_per_worker_unique_results = numel(enc_win_IDs) * numel(maint_win_IDs) * size(unique_pairs, 1) * 8 / (1024^3);  % In GB
        % memory_per_worker_control_EMS_matrices = numel(enc_win_IDs) * numel(maint_win_IDs) * num_test_trials * numSamples * 8 / (1024^3);  % In GB
        % 
        % % Get the number of parallel workers
        % numWorkers = gcp().NumWorkers;
        % 
        % % Estimate total memory usage for all workers
        % totalMemoryUsed = 1 * (memory_per_worker_unique_results + memory_per_worker_control_EMS_matrices);
        % disp(['Estimated total memory usage per worker (GB): ', num2str(totalMemoryUsed)]);

        % %% Fill in control_EMS_matrices using the computed unique results
        % control_EMS_matrices = nan(length(enc_win_IDs), length(maint_win_IDs), num_test_trials, numSamples);
        % % sum([whos().bytes]) / (1024^3) = 5.4875 GB
        % parfor samp = 1:numSamples%parfor samp = 1:numSamples
        %     for j = 1:num_test_trials
        %         % Find the linear index corresponding to (j, samp) in the original pairs_list
        %         linear_idx = (j - 1) * numSamples + samp;
        %         unique_idx = idx_map(linear_idx); % Get the unique pair index from idx_map
        % 
        %         % Assign the precomputed result to the control_EMS_matrices
        %         control_EMS_matrices(:, :, j, samp) = unique_results(:, :, unique_idx);
        %     end
        % end
        % % sum([whos().bytes]) / (1024^3) = 5.4875 GB
        % control_EMS_matrices_std = std(control_EMS_matrices, 0, 3);
        % % control_EMS_matrices = mean(control_EMS_matrices, 3);
        % save(save_file, "chan_id", "control_EMS_matrices", "control_EMS_matrices_std", "chan_test_rows", "-v7.3");
        % all_control_EMS_matrices_by_chan{chan_idx} = control_EMS_matrices;
        % testing memory usage
        % for something = 1:length(all_control_EMS_matrices_by_chan)
        %     all_control_EMS_matrices_by_chan{something} = control_EMS_matrices;
        % end
        clear control_EMS_matrices
        %% legacy: old version for 100 samples without avoiding redudant computations
        % control_EMS_matrices = zeros(length(enc_win_IDs), length(maint_win_IDs), num_trials, numSamples);
        % parfor samp = 1:numSamples
        %     temp_control_EMS_matrices = zeros(length(enc_win_IDs), length(maint_win_IDs), num_trials);
        %     for j = 1:num_test_trials
        %         randtrial = rand_trials(j, samp); % Get the precomputed random trial index
        %         % control_EMS_matrices(:, :, j, samp) = compute_BT_similarity_matrix(chan_test_PSVs(:, :, j), chan_ctrl_PSVs(:, :, randtrial), enc_win_IDs, maint_win_IDs);
        %         temp_control_EMS_matrices(:, :, j) = compute_BT_similarity_matrix(chan_test_PSVs(:, :, j), chan_ctrl_PSVs(:, :, randtrial), enc_win_IDs, maint_win_IDs);
        %     end
        %     control_EMS_matrices(:, :, :, samp) = temp_control_EMS_matrices;
        % end
        % 
        % % Store results in cell arrays for the current channel
        % % all_test_EMS_matrices_by_chan{chan_idx} = test_EMS_matrices;
        % all_control_EMS_matrices_by_chan{chan_idx} = control_EMS_matrices;

        % Save individual channel data
        % save(save_file, "chan_id", "chan_test_table", "chan_ctrl_table", "control_EMS_matrices", "-v7.3");
    end

    %% Save all channel data in a final file
    all_channels_save_file = fullfile(save_folder, 'all_channels_data.mat');
    % save(all_channels_save_file, "item_cor_table", "test_rows", "chan_ctrl_table", "control_table", "all_test_EMS_matrices_by_chan", "all_control_EMS_matrices_by_chan", "unique_channel_IDs", "-v7.3");
    % save(all_channels_save_file, "target_image_ids", "target_enc_ids", "test_rows", "control_rows", "all_control_EMS_matrices_by_chan", "unique_channel_IDs", "-v7.3")
    save(all_channels_save_file, "target_image_ids", "target_enc_ids", "test_rows", "control_rows", "unique_channel_IDs")
end

write_current_script_to_destination(fullfile(params.output_folder, num2str(patient_ID)), mfilename('fullpath'));

% function [label_table, all_window_mean_PS_vectors] = get_cleaned_data(data_file)
%     label_table = data_file.label_table;
%     rows_with_nan = any(isnan(label_table.EMS_means), 2);
%     rows_without_nan = ~rows_with_nan;
%     label_table = label_table(~rows_with_nan,:);
%     all_windowed_mean_PS_vectors = data_file.all_windowed_mean_PS_vectors(:,:,:);
%     all_windowed_mean_PS_vectors = all_windowed_mean_PS_vectors(:,:,rows_without_nan);
% 
%     % remove rows from label_table and all_windowed_mean_PS_vectors if
%     % there are any nans in that row of all_windowed_mean_PS_vectors.
%     reshaped_matrix = reshape(all_windowed_mean_PS_vectors, [], size(all_windowed_mean_PS_vectors, 3));% Reshape the matrix to combine the first two dimensions (891x40) into one
% 
%     nan_vector = any(isnan(reshaped_matrix), 1);  %Check for NaNs along the new combined dimension
%     clear reshaped_matrix
%     nan_vector = nan_vector(:); % Convert the result to a column vector
% 
%     all_window_mean_PS_vectors = all_windowed_mean_PS_vectors(:,:,~nan_vector);
%     label_table = label_table(~nan_vector,:);
% end

function [label_table, all_windowed_mean_PS_vectors] = get_cleaned_data(PS_file)
label_table = PS_file.label_table;
% rows_without_nan = get_rows_without_nan(label_table);

% label_table = label_table(:,:);
all_windowed_mean_PS_vectors = PS_file.all_windowed_mean_PS_vectors(:,:,:);
% all_windowed_mean_PS_vectors = all_windowed_mean_PS_vectors(:,:,rows_without_nan);

end

function rows = filter_rows(label_table, target_image_ids, target_enc_ids)
    % Function to filter rows based on target image and encoding IDs
    % Inputs:
    % - label_table: the table containing the fields 'encoding_correctness' and 'image_id'
    % - target_image_ids: a vector of image IDs to check (e.g., 1:9, [2, 3], etc.)
    % - target_enc_ids: a vector of encoding IDs to check (e.g., 1:3, [1, 3], etc.)

    % Create a logical mask for encoding correctness for the selected target_enc_ids
    enc_correct_mask = ismember(label_table.encoding_correctness(:, target_enc_ids), 1);
    
    % Create a logical mask for matching image IDs in the selected target_enc_ids
    image_id_mask = ismember(label_table.encID_to_imageID(:, target_enc_ids), target_image_ids);

    % Combine both masks using logical AND across the selected target_enc_ids
    rows = any(enc_correct_mask & image_id_mask, 2);
end

function rows = filter_all_correct_without_target_image(label_table, target_image_id)
    % Function to find rows with all encoding IDs correct (1) but excluding target image ID
    % Inputs:
    % - label_table: the table containing the fields 'encoding_correctness' and 'image_id'
    % - target_image_id: the image ID to exclude from all encoding positions

    % Step 1: Create a logical mask for rows where all encoding correctness values are 1
    all_correct_mask = all(label_table.encoding_correctness == 1, 2);

    % Step 2: Create a logical mask for rows where none of the encoding positions contain the target_image_id
    no_target_image_mask = ~any(label_table.encID_to_imageID == target_image_id, 2);

    % Step 3: Combine both masks using logical AND
    rows = all_correct_mask & no_target_image_mask;
end

function rows = filter_by_channel_id(table, target_channel_ID)
    % Function to filter rows based on target channel ID
    % Inputs:
    % - table: the input table containing the field 'channel_ID'
    % - target_channel_ID: the channel ID to filter by

    % Create a logical mask to identify rows with the target channel ID
    rows = table.channel_ID == target_channel_ID;
end

function [filtered_table, filtered_PSVs] = filter_table_PSVs(table, PSVs, rows)
    filtered_PSVs = PSVs(:,:,rows);
    filtered_table = table(rows,:);
end


% function [] = compute_all_between_trial_similarities(data_file)
% 
%     [~, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, ~, non_selection_win_IDs, ~] = get_window_IDs(); % won't change each function call
% 
%     % data_file.all_windowed_mean_PS_vectors: size = all_window_IDs, frequencies, channels x trials 
%     % label_table = data_file.label_table;
% 
%     % whole_trial_ES_matrix = compute_similarity_matrix(mean_PS_vectors, enc1_win_IDs, all_win_IDs);
%     fprintf('\n')
% 
%     %% for each unique channel-item-encID pool the datasubset, remove same-item-wise data from all3enc data, compute
%     % random pairs, process
%     [label_table, all_window_mean_PS_vectors] = get_cleaned_data(data_file);
% 
%     %% pool the target channel's data by channel ID
% 
%     %% pool the target_image_id correct data from channel's data
%     target_image_ids = 1;
%     target_enc_ids = 1;
%     rows = filter_rows(label_table, target_image_ids, target_enc_ids);
% 
%     item_cor_table = label_table(rows,:);
%     item_cor_PSVs = all_window_mean_PS_vectors(:,:,rows);
% 
%     % choose enc window IDs from target_enc_ids
%     switch target_enc_ids
%         case 1
%             enc_win_IDs = enc1_win_IDs;
%         case 2
%             enc_win_IDs = enc2_win_IDs;
%         case 3
%             enc_win_IDs = enc3_win_IDs;
%     end
%     %% generate all3enc correct data without the target_image_id correct data from all3cor_rows
%     % pool the all 3 encoding correct data from target channel's data
%     % rows = label_table.enc_correctness(:,target_encID) == 1; % rows where 3rd item was correct
%     control_rows = filter_all_correct_without_target_image(label_table, target_image_ids);
%     control_table = label_table(control_rows,:);
%     control_PSVs = all_window_mean_PS_vectors(:,:,control_rows);
% 
%     %% compute between-trial random trial pairs between item-present-and-correct encoding trials and correct-but-without-item maintenance of the same channel
%     unique_channel_IDs = unique();
%     save_folder  = fullfile(strrep(data_file.Properties.Source, '.mat', sprintf('all3corrWOitem1_and_item1correctW.mat')));
%     if ~exist(save_folder, 'dir')
%         mkdir(save_folder)
%     end
% 
%     for chan_idx = 1:length(unique_channel_IDs)
%         chan_id = unique_channel_IDs(chan_idx);
%         save_file = fullfile(save_folder, sprintf('BT_%d.mat', chan_id));
%         chan_test_rows = filter_by_channel_id(item_cor_table, chan_id);
%         [chan_test_table, chan_test_PSVs] = filter_table_PSVs(item_cor_table,item_cor_PSVs, chan_test_rows);
% 
%         chan_ctrl_rows = filter_by_channel_id(control_table, chan_id);
%         [chan_ctrl_table, chan_ctrl_PSVs] = filter_table_PSVs(control_table,control_PSVs, chan_ctrl_rows);
% 
%         test_EMS_matrices = zeros(length(enc_win_IDs), length(maint_win_IDs), size(chan_test_PSVs, 3));
%         for j = 1:size(chan_test_PSVS, 3)
%         % compute test matrix
%             test_EMS_matrices(:,:,j) = compute_similarity_matrix(chan_test_PSVs(:,:,j), enc_win_IDs, maint_win_IDs);
%         end
% 
%         % compute control matrices
%         numSamples = 100;
%         control_EMS_matrices = zeros(length(enc_win_IDs), length(maint_win_IDs), size(chan_test_PSVs, 3), numSamples);
%         for samp = 1:numSamples
%             for j = 1:size(chan_test_PSVs, 3)
%                 randtrial = randi(size(chan_ctrl_PSVS, 3));
%                 control_EMS_matrices(:,:,j,samp) = compute_BT_similarity_matrix(chan_test_PSVs(:,:,j), chan_ctrl_PSVs(:,:,randtrial), enc_win_IDs, maint_win_IDs);
%             end
%         end
% 
%     end
%     save(all_channels_save_file, "item_cor_table", "control_table", "test_EMS_matrices_by_chan", "all_control_EMS_matrices_by_chan","unique_channel_IDs")
% end

function pairs_list = filter_trial_pairs_for_within_item(pairs_list, num_ctrl_trials)
    % Filter invalid pairs (test_trial == control_trial) for within_item
    % Find invalid pairs
    invalid_idx = pairs_list(:,1) == pairs_list(:,2);
    
    % Regenerate control_trial values for invalid pairs
    i=0;
    while any(invalid_idx)
        i=i+1;
        if i ==1000
            error('too many attempts to match this trial with another trial. Is it the only available trial?')
        end
        % Replace invalid control_trial values with new random values
        pairs_list(invalid_idx, 2) = randi(num_ctrl_trials, sum(invalid_idx), 1);
        
        % Recompute invalid indices
        invalid_idx = pairs_list(:,1) == pairs_list(:,2);
    end
end