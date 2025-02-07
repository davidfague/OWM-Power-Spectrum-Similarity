%% Now create example of mix-matching trial similarity (i.e. random pairing of trial1 enc with trial2 maint)
% Note that each data subfolder is a single trial's data so we must gather
% trial data. To pool common channel-item-etc data across trials we can
% compare the subfolder names with trialID removed (non_trial_name)
cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
addpath 'Raw Data Storage'               
addpath 'subfunctions'
%% TODO include other encodings.
% Define the folder and the criteria
folder_to_analyze = 'D:\Power Spectrum Similarity\output\201907 gammamod Amyg P046 enc3 correct';
% 'C:\Users\david\Downloads\Power Spectrum Similarity\output\201907 gammamod Amyg P046 enc3 correct';
saving_name = 'between-trial analysis';
folder_to_save_to = strrep(folder_to_analyze, 'output','output1');
folder_to_save_to = fullfile(folder_to_save_to, saving_name);

patients_to_analyze = [201907];
images_to_analyze = {'P046'};
brain_locations_to_analyze = {'amyg'};
channels_to_analyze = {};%{'76'};
trials_to_analyze = {};  % Empty means no filtering on this
encID_to_analyze = 1:3;
encCorrectness_to_analyze = 1;
gamma_chans_to_analyze = [true];

%% window times
time = 1:9001;
% 100ms time window with step of 10ms
wt = length(time);
num_windows = floor((wt - 100) / 10) + 1; % Correct calculation for number of windows
% Precompute window start and end times
window_start_times = time(1:10:(num_windows - 1) * 10 + 1); % Start times
window_end_times = time(window_start_times + 99);            % End times

% access window start and end: (:, window_id)
% windows for 2nd iteration (whole trial windows)
% window_IDs2 = true(size(time));
% window_IDs2 = find_valid_window_IDs_from_ntimes_logical_array(window_IDs2, window_start_times, window_end_times);
% % Define the time windows for each encoding period
% stim1_start = 1000; stim1_end = 1500; % stim 1 time window
% stim2_start = 1500; stim2_end = 2000; % stim 2 time window
% stim3_start = 2000; stim3_end = 2500; % stim 3 time window

%% new windows to compare against
fixation_times = (time >= 0 & time < 1000);
fixation_window_IDs = find_valid_window_IDs_from_ntimes_logical_array(fixation_times, window_start_times, window_end_times);

maintenance_times = (time >= 2500 & time < 6500);
maintenance_window_IDs = find_valid_window_IDs_from_ntimes_logical_array(maintenance_times, window_start_times, window_end_times);

clear time num_windows wt window_start_times window_end_times
%% filter subfolders of folder_to_analyze by criteria
fileList = dir(folder_to_analyze);

% Filter out non-directories and exclude '.' and '..'
dirFlags = [fileList.isdir];  % Get a logical array for directories
subFolders = fileList(dirFlags);  % Keep only the directories
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Exclude '.' and '..'

% Filter the list of folders
subFolders = filter_subfolders(subFolders, patients_to_analyze, brain_locations_to_analyze, channels_to_analyze, gamma_chans_to_analyze, images_to_analyze, trials_to_analyze, encID_to_analyze, encCorrectness_to_analyze);
clear fileList dirFlags patients_to_analyze brain_locations_to_analyze channels_to_analyze gamma_chans_to_analyze images_to_analyze trials_to_analyze encID_to_analyze encCorrectness_to_analyze
%% Get a list non-trial-unique file names and their corresponding files
% Initialize containers for non-trial names and the corresponding original folder names
non_trial_names = {};
original_folders_by_non_trial = containers.Map('KeyType', 'char', 'ValueType', 'any');  % Using a Map to store folder names by unique non-trial name

% Loop over each folder
for k = 1:length(subFolders)
    name = subFolders(k).name;
    if strcmp(name, saving_name)
        continue % have to skip this subfolder
    end
    parts = strsplit(name, '_');  % Split the folder name into parts
    
    % Remove the trial part (assuming trial is always the 6th part)
    non_trial_name = strjoin(parts([1:5, 7:end]), '_');  % Rebuild the name without the trial part (skip part 6)
    
    % If the non-trial name is not already stored, add it to the list
    if ~isKey(original_folders_by_non_trial, non_trial_name)
        non_trial_names{end+1} = non_trial_name;  %#ok<AGROW>
        original_folders_by_non_trial(non_trial_name) = {subFolders(k).name};  % Initialize with the current folder name
    else
        % Retrieve the current list of folders for this non-trial name
        current_folders = original_folders_by_non_trial(non_trial_name);
        
        % Append the new folder name
        current_folders{end+1} = subFolders(k).name;
        
        % Update the map with the modified list
        original_folders_by_non_trial(non_trial_name) = current_folders;
    end
end

clear parts current_folders non_trial_name subFolders

%% analyze between trial data - use non_trial_name pool of trial_folders
data_file_name = 'wholeTrial_ES.mat'; %'wholeTrial_ES_vectors.mat'
for non_trial_name_id = 1:length(non_trial_names) % iterate through channel data
    non_trial_name = non_trial_names{non_trial_name_id}; % common data name without trial indicated
    trial_folders = original_folders_by_non_trial(non_trial_name); % trial data for
    
    if length(trial_folders) < 2
        warning(strcat(non_trial_name, ' has ', num2str(length(trial_folders)), ' trial folders'))
        continue
    end

    % pool data from trial_folders
    all_fixation_mean_PS_vectors = cell(length(trial_folders));
    all_encoding_mean_PS_vectors = cell(length(trial_folders));
    all_maintenance_mean_PS_vectors = cell(length(trial_folders));
    for trial_folder_idx = 1:length(trial_folders)% iterate through trial data
        trial_folder_name = trial_folders{trial_folder_idx};
        % add data from trial_folder to pool
        trial_data_path = fullfile(folder_to_analyze, trial_folder_name, data_file_name);
        trial_data_file = matfile(trial_data_path);

        % get trial data
        encoding_mean_PS_vectors = trial_data_file.window1_mean_PS_vectors;
        fixation_mean_PS_vectors = trial_data_file.window2_mean_PS_vectors(fixation_window_IDs, :);
        maintenance_mean_PS_vectors = trial_data_file.window2_mean_PS_vectors(maintenance_window_IDs, :);

        % add trial data to pool
        all_encoding_mean_PS_vectors{trial_folder_idx} = encoding_mean_PS_vectors;
        all_fixation_mean_PS_vectors{trial_folder_idx} = fixation_mean_PS_vectors;
        all_maintenance_mean_PS_vectors{trial_folder_idx} = maintenance_mean_PS_vectors;

    end
    % save between trial data (save pool of trial data for this
    % non-trial set of data)
    non_trial_save_folder = fullfile(folder_to_save_to, non_trial_name);
    if ~exist(non_trial_save_folder, 'dir')
        mkdir(non_trial_save_folder)
    end
    save_file = fullfile(non_trial_save_folder, 'between_trial_data.mat');
    save(save_file, 'all_maintenance_mean_PS_vectors', 'all_fixation_mean_PS_vectors', 'all_encoding_mean_PS_vectors', 'trial_folders', 'non_trial_name')

    %% compute between-trial EMS - 100 random trial pairs
    num_trial_folders = length(trial_folders);

    if num_trial_folders > 10 % do 100 random pairs
        npairs = 100;
        pairs = zeros(npairs, 2);
        for pair_idx = 1:npairs
            trial_idx1 = randi([1, num_trial_folders]);
            trial_idx2 = trial_idx1; % Initialize num2
            while trial_idx1 == trial_idx2 % get a different trial_idx2
                trial2 = randi([1, randi([1, num_trial_folders]);]);
            end
            % gather the pairs
            pairs(pair_idx,:) = [trial_idx1, trial_idx2];
        end
    else
        npairs = num_trial_folders * num_trial_folders; % Number of all combination
        pairs = zeros(npairs, 2);
        % gather all pairs (less than 100 total)
        pair_idx = 1;
        for trial_idx1 = 1:num_trial_folders
            for trial_idx2 = 1:num_trial_folders
                % Include (i, i) as well
                pairs(pair_idx, :) = [trial_idx1, trial_idx2];
                pair_idx = pair_idx + 1;
            end
        end
    end

    %% do computation across pairs.
    all_trial_data_for_non_trial_name = matfile(savefile);
    all_betweenTrial_EMS_matrices = zeros(length(all_trial_data_for_non_trial_name.all_encoding_mean_PS_vectors{1}), length(all_trial_data_for_non_trial_name.all_maintenance_mean_PS_vectors{1}), size(pairs,1));
    for trial_pair_idx = 1:size(pairs,1)
        enc_trial_idx = pairs(trial_pair_idx, 1); % here idx is the index in trial_folders, which is a subset for this current non_trial_name
        maint_trial_idx = pairs(trial_pair_idx, 2);

        encoding_mean_PS_vectors = all_trial_data_for_non_trial_name.all_encoding_mean_PS_vectors{enc_trial_idx};
        maintenance_mean_PS_vectors = all_trial_data_for_non_trial_name.all_maintenance_mean_PS_vectors{maint_trial_idx};

        betweenTrial_EMS_matrix = compute_similarity_matrix_from_power_vecs(encoding_mean_PS_vectors, maintenance_mean_PS_vectors);

        all_betweenTrial_EMS_matrices(:,:, trial_pair_idx) = betweenTrial_EMS_matrix;
        
    end

end

function similarity_matrix = compute_similarity_matrix_from_power_vecs(mean_power_spectrum_vectors1, mean_power_spectrum_vectors2)
    similarity_matrix = zeros(length(mean_power_spectrum_vectors1), length(mean_power_spectrum_vectors2));
    for i = 1:length(mean_power_spectrum_vectors1)
        window1_mean_power_spectrum_vector = mean_power_spectrum_vectors1(i);
        for j = 1:length(mean_power_spectrum_vectors2)
            window2_mean_power_spectrum_vector = mean_power_spectrum_vectors2(j);
            % Compute correlation and Fisher z-transform
            r = corr(window1_mean_power_spectrum_vector, window2_mean_power_spectrum_vector, 'type', 'spearman');
            z = 0.5 * log((1 + r) / (1 - r));
            % Store similarity value
            similarity_matrix(i, j) = z;
        end
    end
end