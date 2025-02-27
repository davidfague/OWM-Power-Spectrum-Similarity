%% this script is being used to develop the next version of compute_es_btwn_trials.m

% this script will take the output of that script (files for each channel, each file containing a matrix of between-trials similarities.)
% and create a single 'channelMap' using matlab's containerMap, which acts
% as a dictionary from channel to matrix, since the matrices have different
% sizes and cannot fit in one larger matrix.
%% important notes

% reduces storage usage by ~20.75%
% for WI: Total estimated memory used by the map: 436.82 Mb
% for BI: Total estimated memory used by the map: 3816.33 Mb

%% parameters

custom_params.hellbender = false;
custom_params.output_folder_name = 'middle_fixation_baseline';
params = get_parameters(custom_params);
patient_id = 201901;
plot_params.enc_id = 1;
image_id = 1;
params.session_id = 1;

%% example containers.Map usage

% % Create a new empty map
% mapObj = containers.Map('KeyType', 'int32', 'ValueType', 'any');
% 
% % Assign matrices to keys
% mapObj(1) = matrix1;
% mapObj(5) = matrix2;
% mapObj(10) = matrix3;
% 
% % Retrieve a matrix using its key
% retrievedMatrix = mapObj(5);

%% compute map from channel files 
% find the channel files that are within the directory and store their matrices in a containers.Map

% initialize dictionary (containers.Map)
channels_to_btwn_trials_es = containers.Map('KeyType', 'int32' , 'ValueType', 'any');

% get channel_ids_to_use from the channel files location
channel_ids_to_use = extractIntegersFromFilenames(patient_id, params.comp_options, plot_params.enc_id, image_id, params); % image_id=1 is used as a default here. Channels should change with respect to iamge_id.

% get the channel files location (already gotten in
% extractIntegersFromFilenames, but not outputted)
channels_directory = fullfile(sprintf('%s\\%s\\session%d\\%s\\%s\\%d-%d\\enc%s_image%s\\', ...
    params.output_folder, num2str(patient_id), params.session_id, params.comp_options{1}, params.btwn_trial_type,  params.freq_min, params.freq_max, num2str(plot_params.enc_id), num2str(image_id)));
if params.hellbender
    channels_directory = strrep(channels_directory, '\', '/'); % linux instead of windows. I thought fullfile() was supposed to automatically handle that but apparently not.
end

% get similarities for each channel and populate dictionary
for file_idx = 1:length(channel_ids_to_use)
    channel_id = channel_ids_to_use(file_idx);

    channel_file_data = load(sprintf("%s/BT_%d", channels_directory, channel_id));

    channels_to_btwn_trials_es(int32(channel_id)) = channel_file_data;

end

disp(channels_to_btwn_trials_es)

save("testing_channelmap.mat", "channels_to_btwn_trials_es", "-v7.3")

%% compare storage usage of both:
% display storage used in channels_directory
function totalSize = getDirectorySize(directoryPath)
    fileList = dir(fullfile(directoryPath, '**', '*')); % Get all files recursively
    totalSize = sum([fileList.bytes]); % Sum file sizes
    fprintf('Total storage used by directory "%s": %.2f MB\n', directoryPath, totalSize / (1024^2));
end
size_separate_files = getDirectorySize(channels_directory);

% display storaged used by testing_channelmap.mat
function fileSize = getMatFileSize(filePath)
    fileInfo = dir(filePath);
    if isempty(fileInfo)
        error('File does not exist: %s', filePath);
    end
    fileSize = fileInfo.bytes;
    fprintf('Storage used by "%s": %.2f MB\n', filePath, fileSize / (1024^2));
end
size_1_file = getMatFileSize('testing_channelmap.mat');

fprintf("storage used by current method: %.2f MB\n", size_separate_files/ (1024^2))
fprintf("storage used by new method: %.2f MB\n", size_1_file / (1024^2))
fprintf("storage reduction per image, patient, BI, freq_bands: %.2f MB\n",(size_separate_files - size_1_file) / (1024^2))
fprintf("percent storage reduction for BI: %.2f%s\n", ((size_separate_files - size_1_file) / size_separate_files)*100, "%")

%% output: % note that the size of 'all_channels_data.mat' was negligible (929 bytes)
%% BI:

% Total storage used by directory "..\processed_data\middle_fixation_baseline\201901\session1\corr BT BI\EMS\1-40\enc1_image1\": 3243.14 MB
% Storage used by "testing_channelmap.mat": 2565.76 MB

% storage used by current method: 3243.14 MB
% storage used by new method: 2565.76 MB

% storage reduction per image, patient, BI, freq_bands: 677.39 MB
% percent storage reduction for BI: 20.89% MB

%% WI:
% Total storage used by directory "..\processed_data\middle_fixation_baseline\201901\session1\corr BT WI\EMS\1-40\enc1_image1\": 371.02 MB
% Storage used by "testing_channelmap.mat": 294.32 MB

% storage used by current method: 371.02 MB
% storage used by new method: 294.32 MB

% storage reduction per image, patient, BI, freq_bands: 76.70 MB
% percent storage reduction for BI: 20.67%

%% check memory usage of the map:
total_memory = 0;
keys_list = channels_to_btwn_trials_es.keys; % Get all keys
for i = 1:numel(keys_list)
    value = channels_to_btwn_trials_es(keys_list{i});
    mem_info = whos('value');
    total_memory = total_memory + mem_info.bytes;
end
fprintf('Total estimated memory used: %.2f Mb', total_memory /   (1024^2));

% for WI: Total estimated memory used: 436.82 Mb

% for BI: Total estimated memory used: 3816.33 Mb


%% to access the channel:

channels = channels_to_btwn_trials_es.keys;
for chan_idx = 1:length(channels)
    chan_id = channels{chan_idx};
    channel_data = channels_to_btwn_trials_es(chan_id);
end
    