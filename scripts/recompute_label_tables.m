% function label_table = fix_anat_labels(label_table, labelsAnat)
%     channels = unique(label_table.channel_ID);
%     for chan_idx = 1:length(channels)
%         channel_ID = channels(chan_idx);
%         rows_of_channel = label_table.channel_ID == channel_ID;
% 
%         corrected_anat = labelsAnat(channel_ID);
%         previous_anat = unique(label_table.anatomical_label(rows_of_channel));
% 
%         if previous_anat ~= corrected_anat 
%             label_table.anatomical_label(rows_of_channel) = corrected_anat;
%         end
% 
%     end
% end
% 
% function get_preprocessed_data_folder
% 
% end
% 
% function fix_all_label_tables(processed_data_dir)
%     for file = files
%         file_table = load(file, 'label_table');
% 
%         file_table.label_table
%     end
% end


%% Could just alter Zpower_and_PSVs to update them:

custom_params = struct();
custom_params.k12wm = true;
% custom_params.patient_IDs = [201905];
custom_params.baseline_T_lims = [-0.75, -0.25];
custom_params.hellbender = false;
custom_params.output_folder_name = 'middle_fixation_baseline';

params = get_parameters(custom_params);

use_parallel = false; % not implemented
%% loop through patients

for idx = 1:length(params.patient_IDs)
    patient_ID = params.patient_IDs(idx);
    fprintf("Processing patient %s\n", num2str(patient_ID))

    % load patient's preprocessed data
    % patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_ID)));
    patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, patient_ID);

    for session_idx = 1:length(patient_preprocessed_data_paths)
        params.session_id = session_idx;
        patient_preprocessed_data_path = patient_preprocessed_data_paths{session_idx};
        disp(patient_preprocessed_data_path)

        if params.k12wm
            if params.hellbender
                prefix = strsplit(patient_preprocessed_data_path,'/');
            else
                prefix = strsplit(patient_preprocessed_data_path,'\');
            end
            prefix = prefix{end};
            D_OWM_t_file = load(fullfile(patient_preprocessed_data_path, sprintf("%s_1kft_notch_epoch_outliers_bip_demean.mat", prefix)));
            % reformat data
            D_OWM_t_file.D_OWM_t = nan([size(D_OWM_t_file.ftDemean.trial{1,1}),length(D_OWM_t_file.ftDemean.trial)]);
            for i=1:length(D_OWM_t_file.ftDemean.trial)
                D_OWM_t_file.D_OWM_t(:,:,i) = D_OWM_t_file.ftDemean.trial{1,i};
            end
            clear D_OWM_t_file.ftDemean
            % add D_OWM_t_file.labelsanatbkedit
            labels = load(fullfile(patient_preprocessed_data_path, sprintf("%s_labelsAnat.mat", prefix)));
            D_OWM_t_file.labelsanatbkedit = labels.bipolarAnat;
            clear labels
            % load OWM_trial_info
            OWM_trial_info_file = load(fullfile(patient_preprocessed_data_path, sprintf("%s_OWM_trialinfo.mat", prefix)));
            % load the gamma_file
            gamma_file = load(fullfile(patient_preprocessed_data_path, sprintf("%sgammamodchans.mat", prefix)));
        else
            % D_OWM_t_file = matfile(fullfile(patient_path, "D_OWM_t_bipolar.mat"));
            D_OWM_t_file = load(fullfile(patient_preprocessed_data_path, "D_OWM_t_bipolar.mat"));
            OWM_trial_info_file = matfile(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"));
            % load channel gamma modulation (length: nchannels)
            gamma_file = matfile(fullfile(patient_preprocessed_data_path, "gammachans2sd_alltrials.mat"));
    
        end

    
        % store original IDs for labeling purposes
        % nTimes = size(D_OWM_t_file.D_OWM_t, 2); % if not constant
        num_channels = size(D_OWM_t_file.labelsanatbkedit, 1);
    
        num_trials = size(OWM_trial_info_file.Cond_performance, 1);
        orig_channel_IDs = 1:num_channels;
        orig_trial_IDs = 1:num_trials;
    
        % load image info and correctness
        images = OWM_trial_info_file.C; % length 9
        enc_to_imageID = OWM_trial_info_file.Cond_item; % size: (ntrial x 3 enc)
        % if size(enc_to_imageID, 2) == 1 % should be fixed now
        %     enc_to_imageID = reshape(enc_to_imageID, [num_trials, 3]);
        % end

        enc_correctness = OWM_trial_info_file.Cond_performance; % size: (ntrial x 3 enc)
    
        is_gamma_channels = bool_mask_array(gamma_file.sigchans2, num_channels); % converts ([1,3,5], 8) to [1, 0, 1, 0, 1, 0, 0, 0]
    
        clear gamma_file num_channels
        % clear variables no longer needed
    
        %% generate slice indices for the large data subset
        %% get channels to process filtering by brain location and gamma modulation
        % by gamma
        channels_to_process = subset_channels_by_gamma_modulation(orig_channel_IDs, params.use_gamma_mod_chans, is_gamma_channels);
        clear orig_channel_IDs use_gamma_mod_chans
        
        % by anatomical label
        channel_brain_locations = D_OWM_t_file.labelsanatbkedit.anatmacro1;
       
        channels_to_process = subset_channels_by_brain_location(channels_to_process, params.brain_anatomies_to_process, channel_brain_locations);
        channel_brain_locations_to_process = channel_brain_locations(channels_to_process);
    
        %% get trials to process by image, encoding correctness, and encoding ID
        % I think I will process all trials, channels and then subset in post-processing to avoid reprocessing trials.
        trials_to_process = 1:num_trials;
        % by images
        % [image_ids_to_keep, image_to_trialID_encodingID_encodingCorrectness] = subset_trials_by_image_occurence( ...
        %     original_trial_IDs, images_to_process, enc_to_imageID, images, enc_correctness);
        % 
        % % by target encoding ID and correctness
        % image_to_trialID_encodingID_encodingCorrectness = subset_by_encodingCorrectness_and_ID(...
        %     image_ids_to_keep, image_to_trialID_encodingID_encodingCorrectness, target_correctness, target_encodingIDs);
    
        %% generate a label table
        % Create a full grid for channels x trials
        [trial_idx, chan_idx] = ndgrid(1:length(trials_to_process), 1:length(channels_to_process));
        chan_idx = chan_idx(:); % indices of channels_to_process (not original channel IDs)
        trial_idx = trial_idx(:); % indices of trials_to_process (not original trial IDs)
    
        label_table = table( ...
        repmat(patient_ID, length(trials_to_process)*length(channels_to_process), 1), ...                  % Repeat patient ID
        transpose(channels_to_process(chan_idx)), ...                                  % Repeat channel_ID for each trial
        channel_brain_locations_to_process(chan_idx), ...                   % Repeat anatomical_label for each trial
        transpose(is_gamma_channels(channels_to_process(chan_idx))), ...
        enc_to_imageID(trials_to_process(trial_idx),:), ...                            % Repeat encID_to_image for each trial
        enc_correctness(trials_to_process(trial_idx),:), ...                      % Repeat encoding_correctness for each trial
        transpose(trials_to_process(trial_idx)), ...                                    % Repeat trial_ID across channels
        repmat(session_idx, length(trials_to_process)*length(channels_to_process), 1) ... % Repeat session_ID
        );
        clear is_gamma_channels channel_brain_locations_to_process enc_to_imageID enc_correctness
        % Assign meaningful variable names to the table
        label_table.Properties.VariableNames = {'patient_ID', 'channel_ID', 'anatomical_label', 'channel_is_gamma','encID_to_imageID', 'encoding_correctness', 'trial_ID', 'session_ID'};
    
        save_file = get_PS_file(params, num2str(patient_ID), false);
        write_current_script_to_destination(fullfile(params.output_folder, num2str(patient_ID)), strcat(mfilename('fullpath'), '.m'));

        prev_table = load(save_file,'label_table');
        total = length(label_table.anatomical_label);
        old_labels = string(prev_table.label_table.anatomical_label);
        new_labels = string(label_table.anatomical_label);
        same_rows = strcmp(old_labels, new_labels) | (ismissing(old_labels) & ismissing(new_labels));
        compare_table = table(...
            prev_table.label_table.anatomical_label(~same_rows), ...
            label_table.anatomical_label(~same_rows), ...
            'VariableNames', {'old', 'new'});
        percent_new = (1 - (sum(same_rows)/total)) * 100;
        fprintf("percent_new %d\n", percent_new);

        save(save_file, 'label_table', '-v7.3', '-append')
        clear label_table save_file
    end
end

delete(gcp('nocreate'))