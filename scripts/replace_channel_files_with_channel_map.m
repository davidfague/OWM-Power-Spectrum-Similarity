clear all
close all

%% parameters
es_freq_bands = {[1:40]};

custom_params = struct();
custom_params.hellbender = false;
custom_params.output_folder_name = 'middle_fixation_baseline';
% custom_params.patient_IDs = [201907, 201908, 201910, 201915];

target_enc_ids = 1; 
target_image_ids = 1:9;

for k12wm = [true, false]
    custom_params.k12wm = k12wm;
    for i = 1:length(es_freq_bands)
        custom_params.ES_freq_band = es_freq_bands{i}

        params = get_parameters(custom_params); 
        
        %% main loop
        
        for idx = 1:length(params.patient_IDs)
            patient_ID = params.patient_IDs(idx);
            if params.k12wm
                session_ids = get_available_session_ids(params, patient_ID);
            else
                session_ids = 1;
            end
        
            for session_idx = 1:length(session_ids)
                params.session_id = session_idx;
                PS_file = get_PS_file(params, patient_ID, false);
                load(PS_file)
        
                % filter by frequency
                all_windowed_mean_PS_vectors = all_windowed_mean_PS_vectors(:,params.ES_freq_band,:);
                params.freq_min = min(params.ES_freq_band);
                params.freq_max = max(params.ES_freq_band);
            
                % % % compute correlations
                for btwn_trial_type = {'EMS'}%EES
                    params.btwn_trial_type = btwn_trial_type{1};
                    for target_enc_id = target_enc_ids%:3 % can be 1:3W % stim locations, stim 1 vs stim 2 versus stim 3
                        for target_image_id = target_image_ids%[1,3] % can be 1:9
                            for within_item = [true, false] % false is between items
                                params.within_item = within_item;
        
                                fprintf("computing patient%d, enc%d, image%d %d-%dHz\n", patient_ID, target_enc_id, target_image_id, params.freq_min, params.freq_max)
        
                                update_channel_files(patient_ID, target_enc_id, target_image_id, params);
                            end
                        end
                    end
                end
            
                write_current_script_to_destination(fullfile(params.output_folder, num2str(patient_ID)), strcat(mfilename('fullpath'), '.m'));
            
            end
        end
    end
end
%% main computations

function update_channel_files(patient_ID, target_enc_ids, target_image_ids, params)

    % initialize dictionary (containers.Map)
    channels_to_bt_es = containers.Map('KeyType', 'int32' , 'ValueType', 'any');

    channel_ids_to_use = extractIntegersFromFilenames(patient_ID, params.comp_options, target_enc_ids, target_image_ids, params); % image_id=1 is used as a default here. Channels should change with respect to iamge_id.

    if isempty(channel_ids_to_use)
        fprintf(" unable to replace because channel_ids_to_use is empty. RETURNING\n")
        return
    end

    channels_directory = fullfile(sprintf('%s\\%s\\session%d\\%s\\%s\\%d-%d\\enc%s_image%s\\', ...
    params.output_folder, num2str(patient_ID), params.session_id, params.comp_options{1}, params.btwn_trial_type,  params.freq_min, params.freq_max, num2str(target_enc_ids), num2str(target_image_ids)));
    if params.hellbender
        channels_directory = strrep(channels_directory, '\', '/'); % linux instead of windows. I thought fullfile() was supposed to automatically handle that but apparently not.
    end

    save_file = fullfile(channels_directory, sprintf('all_channels_bt_es.mat'));

    % Save all channel data in a final file
    % all_channels_save_file = fullfile(save_folder, 'all_channels_data.mat');
    % save(all_channels_save_file, "item_cor_table", "test_rows", "chan_ctrl_table", "control_table", "all_test_EMS_matrices_by_chan", "all_control_EMS_matrices_by_chan", "unique_channel_IDs", "-v7.3");
    % save(all_channels_save_file, "target_image_ids", "target_enc_ids", "test_rows", "control_rows", "all_control_EMS_matrices_by_chan", "unique_channel_IDs", "-v7.3")
    % save(all_channels_save_file, "target_image_ids", "target_enc_ids", "test_rows", "control_rows", "unique_channel_IDs")

    % get similarities for each channel and populate dictionary
    for file_idx = 1:length(channel_ids_to_use)
        channel_id = channel_ids_to_use(file_idx);
    
        channel_file_data = load(sprintf("%s/BT_%d.mat", channels_directory, channel_id));
    
        channels_to_bt_es(int32(channel_id)) = channel_file_data;

        bt_es = struct();
        bt_es.matrix = channel_file_data.BT_ES;
        bt_es.chan_test_table = channel_file_data.chan_test_table;
        bt_es.chan_ctrl_table = channel_file_data.chan_ctrl_table;

        channels_to_bt_es(int32(channel_id)) = bt_es;

    
    end

    save(save_file, "channels_to_bt_es"); % leaving off here. 
    % Need to check channels_to_bt_es (seems fine), and then
    % delete old
    for file_idx = 1:length(channel_ids_to_use)
        channel_id = channel_ids_to_use(file_idx);
        file_to_delete = fullfile(sprintf('%s/BT_%d.mat', channels_directory, channel_id));
        delete(file_to_delete);
        fprintf(" deleted %s\n", file_to_delete)
    end
    % then update loading this in main_plotting_new

    % channels_to_es_matrix = containers.Map('KeyType', 'int32' , 'ValueType', 'any');

        % bt_es = struct();
        % bt_es_matrix = nan(length(windows1), length(windows2), num_ctrl_trials, num_test_trials); % Initialize results array
        % 
        % % Compute similarity matrices for each test-control trial pair
        % for test_trial_idx = 1:num_test_trials
        %     for ctrl_trial_idx = 1:num_ctrl_trials
        %         if any(any(isnan(chan_test_PSVs(:, :, test_trial_idx)))) || ...
        %            any(any(isnan(chan_ctrl_PSVs(:, :, ctrl_trial_idx)))) || ...
        %            all(all(chan_ctrl_PSVs(:, :, ctrl_trial_idx) == 0)) || ...
        %            all(all(chan_test_PSVs(:, :, test_trial_idx) == 0))
        %             % Skip invalid data
        %             error("skipping because nans or all zeros") % shouldn't get here anymore since above checks. so I replaced with error instead of warning.
        %             bt_es_matrix(:, :, ctrl_trial_idx, test_trial_idx) = nan(size(bt_es_matrix, [1, 2]));
        %         else
        %             % Compute and store the result for the pair
        %             bt_es_matrix(:, :, ctrl_trial_idx, test_trial_idx) = compute_BT_similarity_matrix( ...
        %                 chan_test_PSVs(:, :, test_trial_idx), chan_ctrl_PSVs(:, :, ctrl_trial_idx), ...
        %                 windows1, windows2);
        %         end
        %     end
        % end
        % bt_es.matrix = bt_es_matrix;
        % bt_es.chan_test_table = chan_test_table;
        % bt_es.chan_ctrl_table = chan_ctrl_table;
        % channels_to_es_matrix(int32(chan_id)) = bt_es;
        % 
        % % Save results
        % save(save_file, "bt_es");
        
end