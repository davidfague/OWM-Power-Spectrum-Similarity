%% compute control: correct-item-encoding matched with random non-item, all3 correct maintenance
% note: 'test trials' refer to trials that had encoding taken; 'control trials' refers to trials that had maintenance taken

%% v5 get rid of random pairing & repeated sampling. compute all pairs
%% v4 saving all channels at once? less memory usage may allow since only unique results now
%% v3: use get_parameters.m
%% v2: remove the without_nan subsetting of label_table and PSVs to keep indexing inconsistent; also begin saving only unique results
% also only save the unique trial pairs, do not oversample.
params = get_parameters();
    % temp params
params.within_item = false; % computes control as either within item or between items

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
    PS_file = get_PS_file(params.output_folder, patient_ID, false);
    load(PS_file)

    % % % compute correlations
    for target_enc_ids = 1%:3 % can be 1:3W
        for target_image_ids = [1:9]%[1,3] % can be 1:9
            fprintf("computing patient%d, enc%d, image%d \n", patient_ID, target_enc_ids, target_image_ids)
            [all_channels_save_file] = compute_all_between_trial_similarities(patient_ID, target_enc_ids, target_image_ids, params, label_table, all_windowed_mean_PS_vectors);
        end
    end

    write_current_script_to_destination(fullfile(params.output_folder, num2str(patient_ID)), strcat(mfilename('fullpath'), '.m'));

end

%% main computations
function [all_channels_save_file] = compute_all_between_trial_similarities(patient_ID, target_enc_ids, target_image_ids, params, label_table, all_window_mean_PS_vectors)
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
    if params.within_item % the following could use params.output_folder instead of strrep(PS.source, PS_file, new)
        save_folder = fullfile(sprintf('%s/%d/corr BT WI/enc%d_image%d', params.output_folder, patient_ID, target_enc_ids, target_image_ids));
    else
        save_folder = fullfile(sprintf('%s/%d/corr BT BI/enc%d_image%d', params.output_folder, patient_ID, target_enc_ids, target_image_ids));
    end

    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    %% Loop through each unique channel
    for chan_idx = 1:length(unique_channel_IDs)
        chan_id = unique_channel_IDs(chan_idx);
        save_file = fullfile(save_folder, sprintf('BT_%d.mat', chan_id));
        
        % Filter data for the current channel ID
        chan_test_rows = filter_by_channel_id(item_cor_table, chan_id);
        [chan_test_table, ~] = filter_table_PSVs(item_cor_table, item_cor_PSVs, chan_test_rows);
        [~, chan_test_PSVs] = filter_table_PSVs(item_cor_table, item_cor_PSVs, chan_test_rows);

        chan_ctrl_rows = filter_by_channel_id(control_table, chan_id);
        [chan_ctrl_table, ~] = filter_table_PSVs(control_table, control_PSVs, chan_ctrl_rows);
        [~, chan_ctrl_PSVs] = filter_table_PSVs(control_table, control_PSVs, chan_ctrl_rows);

        num_test_trials = size(chan_test_PSVs, 3);
        num_ctrl_trials = size(chan_ctrl_PSVs, 3);

        %% Generate all possible pairs between test and control trials
        
        % Initialize results array
        BT_ES = nan(length(enc_win_IDs), length(params.maint_win_IDs), num_ctrl_trials, num_test_trials);
    
        % Compute similarity matrices for each test-control trial pair
        fprintf("Computing similarity for %d test trials x %d control trials for chan %d\n", num_test_trials, num_ctrl_trials, chan_id);
        for test_trial_idx = 1:num_test_trials
            for ctrl_trial_idx = 1:num_ctrl_trials
                if any(any(isnan(chan_test_PSVs(:, :, test_trial_idx)))) || ...
                   any(any(isnan(chan_ctrl_PSVs(:, :, ctrl_trial_idx)))) || ...
                   all(all(chan_ctrl_PSVs(:, :, ctrl_trial_idx) == 0)) || ...
                   all(all(chan_test_PSVs(:, :, test_trial_idx) == 0))
                    % Skip invalid data
                    BT_ES(:, :, ctrl_trial_idx, test_trial_idx) = nan(size(BT_ES, [1, 2]));
                else
                    % Compute and store the result for the pair
                    BT_ES(:, :, ctrl_trial_idx, test_trial_idx) = compute_BT_similarity_matrix( ...
                        chan_test_PSVs(:, :, test_trial_idx), chan_ctrl_PSVs(:, :, ctrl_trial_idx), ...
                        enc_win_IDs, params.maint_win_IDs);
                end
            end
        end
    
        % Save results
        save(save_file, "chan_id", "BT_ES", "chan_test_table", "chan_ctrl_table");
        
    end

    %% Save all channel data in a final file
    all_channels_save_file = fullfile(save_folder, 'all_channels_data.mat');
    % save(all_channels_save_file, "item_cor_table", "test_rows", "chan_ctrl_table", "control_table", "all_test_EMS_matrices_by_chan", "all_control_EMS_matrices_by_chan", "unique_channel_IDs", "-v7.3");
    % save(all_channels_save_file, "target_image_ids", "target_enc_ids", "test_rows", "control_rows", "all_control_EMS_matrices_by_chan", "unique_channel_IDs", "-v7.3")
    save(all_channels_save_file, "target_image_ids", "target_enc_ids", "test_rows", "control_rows", "unique_channel_IDs")
end

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

% function [label_table, all_windowed_mean_PS_vectors] = get_cleaned_data(PS_file)
% label_table = PS_file.label_table;
% % rows_without_nan = get_rows_without_nan(label_table);
% 
% % label_table = label_table(:,:);
% all_windowed_mean_PS_vectors = PS_file.all_windowed_mean_PS_vectors(:,:,:);
% % all_windowed_mean_PS_vectors = all_windowed_mean_PS_vectors(:,:,rows_without_nan);
% 
% end

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