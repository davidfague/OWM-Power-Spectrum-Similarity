%% compute fixation-baselined, z-scored power vectors
%% v3 re-enabling bootstrap zbaseline across trials. Need to pass all non-nan channel data instead of individual trial-channel data?
% update save name, update data passed
%% v2: changing file inputs to use OneDrive instead of ExternalHD; use original patient path instead of copying preprocessed data files into 'Raw Data Storage'

%% TODO when recomputing
%1. update saving to use 1 folder per patient, 1 file per large variable:
% done: updated and used get_PS_file
%2. update variable all_windowed_mean_PS_vectors to PSVs
%3. move iter to first dimension of all_windowed_mean_PS_vectors (maybe not
%because will need to adjust elsewhere)

%4. update variable all3_ES_matrix name, save it in correct spot; update
%save_file name

%5. update btwn_trials saving
%% constants
output_folder_name = 'allpatients gammamod allregions allitem allenc baseline across trials';
hellbender = false;
use_gamma_mod_chans = [true];%, false]; % [true, false] means all channels
brain_anatomies_to_process = {};

window_size = 200;

        time = 1:9001; % 1:9001
        wt = length(time);
        num_windows = floor((wt - window_size) / 10) + 1; % Correct calculation for number of windows
        % Precompute window start and end times
        window_start_times = time(1:10:(num_windows - 1) * 10 + 1); % Start times
        window_end_times = time(window_start_times + window_size-1);            % End times

num_frequencies = 100;
% doing all instead of the following specifications
% target_encodingIDs = 1;
% target_correctness = 1;
% images_to_process = {}; %P046

%% Setup Paths

addpath('../subfunctions')

if hellbender
    output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/', output_folder_name); %#ok<UNRCH>
    patient_IDs = [201907, 201908, 201903, 201905, 201906, 201901, 201910, 201915];
else
    output_folder = fullfile('..\AA_Processed Data\', output_folder_name);
    patient_IDs = [201901 201910, 201907]; %#ok<NBRAK>
end
clear output_folder_name

%% Open a local parallel pool if none exists 
if isempty(gcp('nocreate'))
    if hellbender
        parpool('local', 64);%#ok<UNRCH> % more workers
    else
        parpool('local', 2); %#ok<UNRCH>
    end
end

%% loop through patients
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);

    if hellbender
        patient_path = sprintf('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/CS%s/', num2str(patient_ID)); %#ok<UNRCH>
    else
        patient_path = sprintf('../../../..//OWM Utah Data/CS%s/', num2str(patient_ID));
    end

    % D_OWM_t_file = matfile(fullfile(patient_path, "D_OWM_t_bipolar.mat"));
    D_OWM_t_file = load(fullfile(patient_path, "D_OWM_t_bipolar.mat"));
    OWM_trial_info_file = matfile(fullfile(patient_path, "OWM_trialinfo.mat"));

    
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
    gamma_file = matfile(fullfile(patient_path, "gammachans2sd_alltrials.mat"));

    is_gamma_channels = bool_mask_array(gamma_file.sigchans2, num_channels); % converts ([1,3,5], 8) to [1, 0, 1, 0, 1, 0, 0, 0]

    clear gamma_file num_channels
    % clear variables no longer needed

    %% generate slice indices for the large data subset
    %% get channels to process filtering by brain location and gamma modulation
    % by gamma
    channels_to_process = subset_channels_by_gamma_modulation(orig_channel_IDs, use_gamma_mod_chans, is_gamma_channels);
    clear orig_channel_IDs use_gamma_mod_chans
    
    % by anatomical label
    channel_brain_locations = D_OWM_t_file.labelsanatbkedit;
    channel_brain_locations = channel_brain_locations.anatmacro1;
    channels_to_process = subset_channels_by_brain_location(channels_to_process, brain_anatomies_to_process, channel_brain_locations);
    channel_brain_locations_to_process = channel_brain_locations(channels_to_process);

    %% get trials to process by image, encoding correctness, and encoding ID
    %% I think I will process all trials, channels and then subset in post-processing to avoid reprocessing entire trials.
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
    %% new (testing)
    size(label_table)
    length(channels_to_process) == length(unique(label_table.channel_ID))

    all_windowed_mean_PS_vectors = nan(881, 100, size(label_table,1)); % num_windows, num_frequencies, num_rows

    for chan_id = 1:length(channels_to_process)
        original_channel_id = channels_to_process(chan_id);
        % compute Zbaseline power for all trials for this channel
        valid_trials_mask = all(~isnan(D_OWM_t_file.D_OWM_t(chan, :, :)), 2);
        valid_trials = find(valid_trials_mask);

        % Extract all valid trials for the channel
        time_data = D_OWM_t_file.D_OWM_t(chan, :, valid_trials); % Pass all valid trials for this channel to Zpower

        % Compute Zpower for all valid trials
        Zpower = compute_Zpower_v2(time_data);

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
            windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(trial_Zpower);
    
            % Store the result in the preallocated matrix
            iter = iter_map(chan, trial_id);
            all_windowed_mean_PS_vectors(:, :, row_id) = windowed_mean_PS_vectors;

        end
    end


    %% this might have caused a mismatch between label_table and the computed PSVs
    % remove rows that have nan
    rows_with_nan = any(isnan(label_table.encoding_correctness), 2); % remove aborted trials from label_table
    label_table = label_table(~rows_with_nan, :);
    chan_idx = chan_idx(~rows_with_nan); % should correspond to label_table.trial_ID
    trial_idx = trial_idx(~rows_with_nan); % should correspond to label_table.trial_ID

    %% baseline across all valid trials for each channel (new)
    tic;
    % Initialize iter mapping
    num_channels = size(D_OWM_t_file.D_OWM_t, 1);
    num_trials = size(D_OWM_t_file.D_OWM_t, 3);
    
    % Precompute a linear mapping for valid channel-trial combinations
    iter_map = zeros(num_channels, num_trials);
    current_iter = 1;
    
    % Precompute valid channel-trial combinations
    valid_indices = [];
    for chan = 1:num_channels
        valid_trials_mask = all(~isnan(D_OWM_t_file.D_OWM_t(chan, :, :)), 2); % Identify valid trials for this channel
        valid_trials = find(valid_trials_mask); % Indices of valid trials
        valid_indices = [valid_indices; [repmat(chan, length(valid_trials), 1), valid_trials]]; %#ok<AGROW>
    end

    for chan = 1:num_channels
        valid_trials_mask = all(~isnan(D_OWM_t_file.D_OWM_t(chan, :, :)), 2);
        valid_trials = find(valid_trials_mask);
    
        % Assign a unique iter index to each valid channel-trial pair
        iter_map(chan, valid_trials) = current_iter:(current_iter + length(valid_trials) - 1);
        current_iter = current_iter + length(valid_trials);
    end
    
    % Total number of valid combinations
    num_valid_combinations = current_iter - 1;
    
    % Preallocate the result matrix
    all_windowed_mean_PS_vectors = zeros(num_windows, num_frequencies, num_valid_combinations);

    % Adjust label_table if required
    % Use logical indexing to filter label_table
    % is_valid_trial = ismember(label_table.trial_ID, valid_indices(:, 2));
    % label_table = label_table(is_valid_trial, :);
    % Use valid_indices to filter label_table
    is_valid_combination = ismember([chan_idx, trial_idx], valid_indices, 'rows');
    label_table = label_table(is_valid_combination, :);

    % parfor loop
    for chan = 1:num_channels
        % Identify valid trials for this channel
        valid_trials_mask = all(~isnan(D_OWM_t_file.D_OWM_t(chan, :, :)), 2);
        valid_trials = find(valid_trials_mask);
    
        if isempty(valid_trials)
            continue; % Skip this channel if there are no valid trials
        end
    
        % Extract all valid trials for the channel
        time_data = D_OWM_t_file.D_OWM_t(chan, :, valid_trials); % Pass all valid trials for this channel to Zpower
    
        % Compute Zpower for all valid trials
        Zpower = compute_Zpower_v2(time_data);
    
        % Compute windowed mean power spectral vectors for each trial
        for trial_i = 1:length(valid_trials)
            trial_id = valid_trials(trial_i);
    
            % Extract Zpower for the specific trial
            trial_Zpower = squeeze(Zpower(trial_i, :, :)); % Slice the Zpower for the trial
    
            % Compute windowed mean power vectors for the trial
            windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(trial_Zpower);
    
            % Store the result in the preallocated matrix
            iter = iter_map(chan, trial_id);
            all_windowed_mean_PS_vectors(:, :, iter) = windowed_mean_PS_vectors;

        end
    end

    % %% Process individual data slice ( old)
    % % get data subset
    % % data_subset = raw_data_file.D_OWM_t(:, :, :);
    % % data_subset = data_subset(channels_to_process, :, trials_to_process);
    % 
    % % channel-trial combinations and calculate whole-trial windowed-mean power
    % % spectrums
    % % Define the number of iterations
    % num_iterations = length(label_table.patient_ID);
    % % Prepare scalar variables for each iteration to avoid broadcasting
    % chan_idx_list = channels_to_process(chan_idx);  % Precompute the channel indices
    % trial_idx_list = trials_to_process(trial_idx);  % Precompute the trial indices
    % 
    % % Initialize a placeholder for the all_windowed_mean_PS_vectors matrix
    % % by calclating the first one
    % tic;
    % all_windowed_mean_PS_vectors = preallocate_all_windowed_mean_PS_vectors(...
    %     D_OWM_t_file.D_OWM_t(channels_to_process(chan_idx_list(1)), :, trials_to_process(trial_idx_list(1))), ...
    %     num_iterations);
    % % whos all_windowed_mean_PS_vectors % should be around 4 GB
    % parfor iter = 2:num_iterations
    %     %% load batch slice of patient data
    %     % Load the individual data slice for this channel and trial
    %     % Get the specific channel and trial ID for the current iteration
    %     chan_id = chan_idx_list(iter);  % Scalar value of channel index
    %     trial_id = trial_idx_list(iter);  % Scalar value of trial index
    %     time_data = D_OWM_t_file.D_OWM_t(chan_id, :, trial_id);  % Efficient access using matfile
    % 
    %     % Compute Zpower from the time data
    %     Zpower = compute_Zpower_v2(time_data);
    % 
    %     % Compute windowed_mean_PS_vectors
    %     windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(Zpower);
    % 
    %     % Store the current result
    %     all_windowed_mean_PS_vectors(:, :, iter) = windowed_mean_PS_vectors;
    % end
    elapsed=toc;
    fprintf("Finished in %.2f seconds\n", elapsed);
    clear time_data Zpower windowed_mean_PS_vectors
    % save
    % save_file = strcat(fullfile(output_folder, num2str(patient_ID)), '.mat');
    if ~exist(fullfile(output_folder, num2str(patient_ID)), "file")
        mkdir(fullfile(output_folder, num2str(patient_ID)))
    end

    save_file = get_PS_file(output_folder, num2str(patient_ID), false);
    write_current_script_to_destination(fullfile(output_folder, num2str(patient_ID), 'Zpower and PSVs.m'))
    save(save_file, 'all_windowed_mean_PS_vectors', 'label_table', '-v7.3')
    clear all_windowed_mean_PS_vectors label_table save_file
    %% when rerunning could save label_table & PSVs like this:
    % save_folder = fullfile(output_folder, num2str(patient_ID)));
    % table_save_file = fullfile(save_folder, 'label_table.mat')
    % save(table_save_file, 'label_table')
    % PSVs_save_file = fullfile(save_folder, 'PSVs.mat')
    % save(PSVs_save_file, 'all_windowed_mean_PS_vectors
end

delete(gcp('nocreate'))

function all_windowed_mean_PS_vectors = preallocate_all_windowed_mean_PS_vectors(time_data, num_iterations)
        % Compute Zpower from the time data
        Zpower = compute_Zpower(time_data);
        
        % Compute windowed_mean_PS_vectors
        windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(Zpower);

        [windows, frequencies] = size(windowed_mean_PS_vectors);
        all_windowed_mean_PS_vectors = zeros(windows, frequencies, num_iterations);
        all_windowed_mean_PS_vectors(:, :, 1) = windowed_mean_PS_vectors;
end

function windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(Zpower)
    % Declare persistent variables so that we can define them only once and
    % avoid passing them.
    persistent time wt num_windows window_start_times window_end_times ES_freq_band window_IDs_to_use window_size
    % Check if persistent variables are empty (first function call)
    if isempty(time)
        % These values are computed once and reused
        time = 1:9001; % 1:9001
        window_size=200;
        wt = length(time);
        num_windows = floor((wt - window_size) / 10) + 1; % Correct calculation for number of windows
        % Precompute window start and end times
        window_start_times = time(1:10:(num_windows - 1) * 10 + 1); % Start times
        window_end_times = time(window_start_times + window_size-1);             % End times
        % Initialize and process the window IDs once
        window_IDs_to_use = true(size(time)); % Choose all times initially
        window_IDs_to_use = find_valid_window_IDs_from_ntimes_logical_array(window_IDs_to_use, window_start_times, window_end_times); % Filter valid IDs
        % Define the frequency band
        ES_freq_band = 1:100;
    end

    % Preallocate the output matrix
    windowed_mean_PS_vectors = zeros(length(window_IDs_to_use), length(ES_freq_band));
    
    % Compute the windowed mean power spectral vectors
    for i = 1:length(window_IDs_to_use)
        window_ID = window_IDs_to_use(i);
        windowed_mean_PS_vectors(i, :) = mean(Zpower(ES_freq_band, window_start_times(window_ID):window_end_times(window_ID)), 2)';
    end
end