%% compute fixation-baselined, z-scored power vectors
%% v5 implementing parameters.m
%% v4 redoing the new v3 computation because getting the 'valid' indices & subsetting was becoming too complicated and not working when trying to pass only valid trials to compute_Zpower.
% now saving PSVs and label_table with nrows = (trials_to_process x channels_to_process) instead
% of doing some extra filtering after initialization and before processing.
%% v3 re-enabling bootstrap zbaseline across trials. Need to pass all non-nan channel data instead of individual trial-channel data?
% update save name, update data passed
%% v2: changing file inputs to use OneDrive instead of ExternalHD; use original patient path instead of copying preprocessed data files into 'Raw Data Storage'
%% Setup Paths
params = get_parameters();
addpath('../subfunctions')

%% Open a local parallel pool if none exists 
if isempty(gcp('nocreate'))
    if params.hellbender
        parpool('local', params.num_wrkrs_hb);
    else
        parpool('local', params.num_wrkrs_local); 
    end
end

%% loop through patients
for idx = 1:length(params.patient_IDs)
    patient_ID = params.patient_IDs(idx);

    % load patient's preprocessed data
    patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_ID)));

    % D_OWM_t_file = matfile(fullfile(patient_path, "D_OWM_t_bipolar.mat"));
    D_OWM_t_file = load(fullfile(patient_preprocessed_data_path, "D_OWM_t_bipolar.mat"));
    OWM_trial_info_file = matfile(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"));

    
    % store original IDs for labeling purposes
    % nTimes = size(D_OWM_t_file.D_OWM_t, 2); % if not constant
    num_channels = size(D_OWM_t_file.labelsanatbkedit, 1);

    num_trials = size(OWM_trial_info_file.Cond_item, 1);
    orig_channel_IDs = 1:num_channels;
    orig_trial_IDs = 1:num_trials;

    % load image info and correctness
    images = OWM_trial_info_file.C; % length 9
    enc_to_imageID = OWM_trial_info_file.Cond_item; % size: (ntrial x 3 enc)
    enc_correctness = OWM_trial_info_file.Cond_performance; % size: (ntrial x 3 enc)

    % load channel gamma modulation (length: nchannels)
    gamma_file = matfile(fullfile(patient_preprocessed_data_path, "gammachans2sd_alltrials.mat"));

    is_gamma_channels = bool_mask_array(gamma_file.sigchans2, num_channels); % converts ([1,3,5], 8) to [1, 0, 1, 0, 1, 0, 0, 0]

    clear gamma_file num_channels
    % clear variables no longer needed

    %% generate slice indices for the large data subset
    %% get channels to process filtering by brain location and gamma modulation
    % by gamma
    channels_to_process = subset_channels_by_gamma_modulation(orig_channel_IDs, params.use_gamma_mod_chans, is_gamma_channels);
    clear orig_channel_IDs use_gamma_mod_chans
    
    % by anatomical label
    channel_brain_locations = D_OWM_t_file.labelsanatbkedit;
    channel_brain_locations = channel_brain_locations.anatmacro1;
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
    transpose(trials_to_process(trial_idx)) ...                                    % Repeat trial_ID across channels
    );
    clear is_gamma_channels channel_brain_locations_to_process enc_to_imageID enc_correctness
    % Assign meaningful variable names to the table
    label_table.Properties.VariableNames = {'patient_ID', 'channel_ID', 'anatomical_label', 'channel_is_gamma','encID_to_imageID', 'encoding_correctness', 'trial_ID'};
    
    %% computing
    tic;
    if ~(length(channels_to_process) == length(unique(label_table.channel_ID)))
        error("channels_to_process should be the unique channels in label_table, [%s %s]", num2str(length(channels_to_process)), numstr(length(unique(label_table.channel_ID))))
    end

    all_windowed_mean_PS_vectors = nan(params.num_windows, params.num_PSV_frequencies, size(label_table,1)); % num_windows, num_frequencies, num_rows

    for chan_id = 1:length(channels_to_process)
        original_channel_id = channels_to_process(chan_id);
        % compute Zbaseline power for all trials for this channel
        valid_trials_mask = all(~isnan(D_OWM_t_file.D_OWM_t(original_channel_id, :, :)), 2);
        valid_trials = find(valid_trials_mask);

        % Extract all valid trials for the channel
        time_data = D_OWM_t_file.D_OWM_t(original_channel_id, :, valid_trials); % Pass all valid trials for this channel to Zpower

        % Compute Zpower for all valid trials
        Zpower = compute_Zpower_v3(time_data, params);

        % Compute windowed mean power spectral vectors for each trial
        for trial_i = 1:length(valid_trials)
            trial_id = valid_trials(trial_i);

            %find row of 
            row_id = find(label_table.channel_ID==original_channel_id & label_table.trial_ID == trial_id);
            if length(row_id) > 1
                error("should be one row per chan x trial. %s pairs found", num2str(length(row_id)))
            end
    
            % Extract Zpower for the specific trial
            trial_Zpower = squeeze(Zpower(trial_i, :, :)); % Slice the Zpower for the trial
    
            % Compute windowed mean power vectors for the trial
            windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(trial_Zpower, params);
    
            % Store the result in the preallocated matrix
            all_windowed_mean_PS_vectors(:, :, row_id) = windowed_mean_PS_vectors;

        end
    end

    elapsed=toc;
    fprintf("Finished in %.2f seconds\n", elapsed);
    clear time_data Zpower windowed_mean_PS_vectors
    % save
    % save_file = strcat(fullfile(output_folder, num2str(patient_ID)), '.mat');
    if ~exist(fullfile(output_folder, num2str(patient_ID)), "file")
        mkdir(fullfile(output_folder, num2str(patient_ID)))
    end

    save_file = get_PS_file(output_folder, num2str(patient_ID), false);
    write_current_script_to_destination(fullfile(output_folder, num2str(patient_ID)), mfilename('fullpath'));
    save(save_file, 'all_windowed_mean_PS_vectors', 'label_table', '-v7.3')
    clear all_windowed_mean_PS_vectors label_table save_file
end

delete(gcp('nocreate'))

function windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(Zpower, params)

    % Preallocate the output matrix
    windowed_mean_PS_vectors = zeros(length(params.window_IDs_to_use), length(params.PSV_freq_band));
    
    % Compute the windowed mean power spectral vectors
    for i = 1:length(params.window_IDs_to_use)
        window_ID = params.window_IDs_to_use(i);
        windowed_mean_PS_vectors(i, :) = mean(Zpower(params.PSV_freq_band, params.window_start_times(window_ID):params.window_end_times(window_ID)), 2)';
    end
end