%% compute fixation-baselined, z-scored power vectors

%% TODO when recomputing
%1. update saving to use 1 folder per patient, 1 file per large variable
%2. update variable all_windowed_mean_PS_vectors to PSVs
%3. move iter to first dimension of all_windowed_mean_PS_vectors (maybe not
%because will need to adjust elsewhere)

%4. update variable all3_ES_matrix name, save it in correct spot; update
%save_file name

%5. update btwn_trials saving
%% constants
hellbender = false;
use_gamma_mod_chans = [true];%, false]; % [true, false] means all channels
brain_anatomies_to_process = {};
% target_encodingIDs = 1;
% target_correctness = 1;
% images_to_process = {}; %P046
%%
% cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
if hellbender
    addpath 'Raw Data Storage'               %#ok<UNRCH>
    addpath 'subfunctions'
else
    addpath('../Z_Raw Data Storage')      %#ok<UNRCH>    
    addpath('../subfunctions')
end

%% Specify large data subset to process - empty means all - regions, gamma/non-gamma channels, trials, performance, items
if hellbender
    patient_IDs = [201907 201908, 201903, 201905, 201906, 201901, 201910, 201915]; %#ok<UNRCH>
else
    patient_IDs = [201907]; %#ok<NBRAK2>
end

% output_folder = 'D:\Power Spectrum Similarity\parallel output\allpatients gammamod allregions allitem enc1 correct';
if hellbender
    output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc'); %#ok<UNRCH>
else 
    % output_folder = 'D:\Power Spectrum Similarity\parallel
    % output\allpatients gammamod allregions allitem enc1 correct'; 
    output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct'; %#ok<UNRCH>
end

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

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

    %% prepare patient's raw data
    % enable partial load of raw data
    if hellbender
        raw_data_file = matfile(fullfile('Raw Data Storage', ['D_OWM_t_bipolar_', num2str(patient_ID), '.mat'])); %#ok<UNRCH>    
        raw_info_file = matfile(fullfile('Raw Data Storage', ['OWM_trialinfo_', num2str(patient_ID), '.mat']));
    else
        raw_data_file = matfile(fullfile('Z_Raw Data Storage', ['D_OWM_t_bipolar_', num2str(patient_ID), '.mat'])); %#ok<UNRCH>    
        raw_info_file = matfile(fullfile('Z_Raw Data Storage', ['OWM_trialinfo_', num2str(patient_ID), '.mat']));
    end
    
    % store original IDs for labeling purposes
    % nTimes = size(raw_data_file.D_OWM_t, 2); % if not constant
    num_channels = size(raw_data_file.labelsanatbkedit, 1);

    num_trials = size(raw_info_file.Cond_item, 1);
    orig_channel_IDs = 1:num_channels;
    orig_trial_IDs = 1:num_trials;

    % load image info and correctness
    images = raw_info_file.C; % length 9
    enc_to_imageID = raw_info_file.Cond_item; % size: (ntrial x 3 enc)
    enc_correctness = raw_info_file.Cond_performance; % size: (ntrial x 3 enc)

    % load channel gamma modulation (length: nchannels)
    if hellbender
        gamma_data = matfile(fullfile('Raw Data Storage', ['gammachans2sd_alltrials_', num2str(patient_ID), '.mat'])); %#ok<UNRCH>
    else
        gamma_data = matfile(fullfile('Z_Raw Data Storage', ['gammachans2sd_alltrials_', num2str(patient_ID), '.mat'])); %#ok<UNRCH>
    end
    is_gamma_channels = bool_mask_array(gamma_data.sigchans2, num_channels); % converts ([1,3,5], 8) to [1, 0, 1, 0, 1, 0, 0, 0]

    clear gamma_data num_channels
    % clear variables no longer needed

    %% generate slice indices for the large data subset
    %% get channels to process filtering by brain location and gamma modulation
    % by gamma
    channels_to_process = subset_channels_by_gamma_modulation(orig_channel_IDs, use_gamma_mod_chans, is_gamma_channels);
    clear orig_channel_IDs use_gamma_mod_chans
    
    % by anatomical label
    channel_brain_locations = raw_data_file.labelsanatbkedit;
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

    %% this might have caused a mismatch between label_table and the computed PSVs
    % remove rows that have nan
    rows_with_nan = any(isnan(label_table.encoding_correctness), 2); % remove aborted trials from label_table
    label_table = label_table(~rows_with_nan, :);
    chan_idx = chan_idx(~rows_with_nan);
    trial_idx = trial_idx(~rows_with_nan);

    %% Process individual data
    % get data subset
    % data_subset = raw_data_file.D_OWM_t(:, :, :);
    % data_subset = data_subset(channels_to_process, :, trials_to_process);

    % channel-trial combinations and calculate whole-trial windowed-mean power
    % spectrums
    % Define the number of iterations
    num_iterations = length(label_table.patient_ID);
    
    % Initialize a placeholder for the all_windowed_mean_PS_vectors matrix
    all_windowed_mean_PS_vectors = preallocate_all_windowed_mean_PS_vectors(...
        raw_data_file.D_OWM_t(channels_to_process(chan_idx_list(1)), :, trials_to_process(trial_idx_list(1))), ...
        num_iterations);

    % whos all_windowed_mean_PS_vectors % should be around 4 GB

    % Prepare scalar variables for each iteration to avoid broadcasting
    chan_idx_list = channels_to_process(chan_idx);  % Precompute the channel indices
    trial_idx_list = trials_to_process(trial_idx);  % Precompute the trial indices
    for iter = 2:num_iterations
        %% load batch slice of patient data
        % Load the individual data slice for this channel and trial
        % Get the specific channel and trial ID for the current iteration
        chan_id = chan_idx_list(iter);  % Scalar value of channel index
        trial_id = trial_idx_list(iter);  % Scalar value of trial index
        time_data = raw_data_file.D_OWM_t(chan_id, :, trial_id);  % Efficient access using matfile
        
        % Compute Zpower from the time data
        Zpower = compute_Zpower(time_data);
        
        % Compute windowed_mean_PS_vectors
        windowed_mean_PS_vectors = compute_windowed_mean_PS_vectors(Zpower);
        
        % Store the current result
        all_windowed_mean_PS_vectors(:, :, iter) = windowed_mean_PS_vectors;
    end
    clear time_data Zpower windowed_mean_PS_vectors
    % save
    save_file = strcat(fullfile(output_folder, num2str(patient_ID)), '.mat');
    save(save_file, 'all_windowed_mean_PS_vectors', 'label_table', '-v7.3')
    clear all_windowed_mean_PS_vectors label_table
    %% when rerunning want to save like this:
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
    persistent time wt num_windows window_start_times window_end_times ES_freq_band window_IDs_to_use
    % Check if persistent variables are empty (first function call)
    if isempty(time)
        % These values are computed once and reused
        time = 1:size(Zpower, 2); % 1:9001
        wt = length(time);
        num_windows = floor((wt - 100) / 10) + 1; % Correct calculation for number of windows
        % Precompute window start and end times
        window_start_times = time(1:10:(num_windows - 1) * 10 + 1); % Start times
        window_end_times = time(window_start_times + 99);            % End times
        % Initialize and process the window IDs once
        window_IDs_to_use = true(size(time)); % Choose all times initially
        window_IDs_to_use = find_valid_window_IDs_from_ntimes_logical_array(window_IDs_to_use, window_start_times, window_end_times); % Filter valid IDs
        % Define the frequency band
        ES_freq_band = 1:40;
    end

    % Preallocate the output matrix
    windowed_mean_PS_vectors = zeros(length(window_IDs_to_use), length(ES_freq_band));
    
    % Compute the windowed mean power spectral vectors
    for i = 1:length(window_IDs_to_use)
        window_ID = window_IDs_to_use(i);
        windowed_mean_PS_vectors(i, :) = mean(Zpower(ES_freq_band, window_start_times(window_ID):window_end_times(window_ID)), 2)';
    end
end