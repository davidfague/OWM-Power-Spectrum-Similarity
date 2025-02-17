% TODO: add date to when saving this script? add more specific location for
% each time running?
% params.m to keep things consistent across files for easy parameter variation

function params = get_parameters(custom_params, save_script, save_values)

%% specify working directory

% check local directory below too
working_dir = 'C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\OWM-Power-Spectrum-Similarity\'; % maybe call it repo_dir instead. Actual working_dir would be /scripts
working_dir = "/mnt/pixstor/data/drfrbc/OWM-Power-Spectrum-Similarity/";
%% addpath all utils

% addpath(fullfile(strcat(working_dir,'\subfunctions')))
utils_dir = fullfile(strcat(working_dir,'utils'));
addSubDirs(utils_dir)

function addSubDirs(dirName)
    % Adds the given directory and all its subdirectories (except those
    % named '.', '..', and '.git') to the MATLAB path.

    if ~isfolder(dirName) % Check that the directory exists.
        error('Directory does not exist: %s', dirName);
    end

    addpath(dirName); % Add the current directory to the MATLAB path.

    items = dir(dirName); % Get a list of all items in the directory.
    for k = 1:length(items) % Loop over each item.
        if items(k).isdir % Check if the item is a directory.
            if ~ismember(items(k).name, {'.', '..', '.git'}) % Exclude '.', '..', and '.git'.
                subDir = fullfile(dirName, items(k).name); % Build the full path for the subdirectory.
                addSubDirs(subDir); % Recursively add subdirectories.
            end
        end
    end
end

%% defaults 

params = struct();
params.k12wm = false;          % Set to true for k12wm-specific settings
params.hellbender = false;     % true for supercomputer; false for local
params.num_wrkrs_hb = 64;
params.num_wrkrs_local = 2;
params.init_par_pool = false;  % enable to run parfor
params.output_folder_name = 'allpatients gammamod allregions allitem allenc baseline across trials';
params.processed_data_dir = '../processed_data/';
params.ES_freq_band = 1:40;

%% Override first set of defaults with any custom parameters

% example usage:
    % custom_params.num_wrkrs_hb = 128;
    % custom_params.output_folder_name = 'modified_output_folder';
    % params = get_parameters(custom_params);

if nargin >= 1 && ~isempty(custom_params)
    fnames = fieldnames(custom_params);
    for i = 1:numel(fnames)
        % % If the custom field exists, override the default value.
        % % (Extra fields not already in params produce a warning.)
        % if isfield(params, fnames{i})
            params.(fnames{i}) = custom_params.(fnames{i});
        % else
        %     warning('Unknown parameter: %s. Ignored.', fnames{i});
        % end
    end
end

%% logic and extras

if params.hellbender
    params.patient_IDs = [201901, 201902, 201903, 201905, 201906, 201907, 201908, 201910, 201915];
    params.output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/pprocessed_data/', params.output_folder_name);
    params.preprocessed_data_location = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/');
    params.local_directory = fullfile('/home/drfrbc/Power Spectrum Similarity/');
else
    params.patient_IDs = [201901, 201902, 201903, 201905, 201906, 201907, 201908, 201910, 201915];
    params.output_folder = fullfile(params.processed_data_dir, params.output_folder_name);
    params.preprocessed_data_location = fullfile('../../../..//OWM Utah Data/');
    params.local_directory = fullfile('C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\OWM-Power-Spectrum-Similarity\');
end

if params.k12wm
    params.time = 1:7001;
    params.patient_IDs = [004, 005, 006, 007, 008, 009, 010];
    if params.hellbender
        params.output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/pprocessed_data/', params.output_folder_name, '/k12wm/');
        params.preprocessed_data_location = fullfile('/cluster/VAST/bkybg-lab/Data/k12wm/');
    else
        params.output_folder = fullfile(params.processed_data_dir, params.output_folder_name, '\k12wm\');
        params.preprocessed_data_location = fullfile('../../../..//k12wm/');
        params.local_directory = fullfile('C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\OWM-Power-Spectrum-Similarity\');
    end
else
    params.time = 1:9001;
end

%% hyper-parameters for processing

% Define the time windows for each encoding period
params.stim1_start = 1000; params.stim1_end = 1500; % stim 1 time window
params.stim2_start = 1500; params.stim2_end = 2000; % stim 2 time window
params.stim3_start = 2000; params.stim3_end = 2500; % stim 3 time window
params.maintenance_duration = 4000;

% windows for computing PSVs
params.window_size = 200;
wt = length(params.time); % time from preprecess data values
params.num_windows = floor((wt - params.window_size) / 10) + 1; % Correct calculation for number of windows
params.window_start_times = params.time(1:10:(params.num_windows - 1) * 10 + 1); % Start times
params.window_end_times = params.time(params.window_start_times + params.window_size-1);             % End times
params.window_IDs_to_use = true(size(params.time)); % Choose all times initially for PSVs
params.window_IDs_to_use = find_valid_window_IDs_from_ntimes_logical_array(params.window_IDs_to_use, params.window_start_times, params.window_end_times); % Filter valid IDs

% windows for computing similarity between PSV subsets
[params.fixation_win_IDs, params.enc1_win_IDs, params.enc2_win_IDs, params.enc3_win_IDs, params.maint_win_IDs, params.non_recall_win_IDs, params.all_win_IDs] = get_window_IDs(params);

% for computing power
params.PSV_freq_band = 1:100;
params.num_PSV_frequencies = length(params.PSV_freq_band);
params.Fsample = 1000;
params.wavenum = 6;

% baselining Zpower
params.baseline_across_trials = true; % false does within trial
params.baseline_T_lims = [-0.75, -0.25];% considering T in seconds and T=T-1 so that 0=stim onset

% data subsetting on initial Zpower computation
params.brain_anatomies_to_process = {}; % only used in first script: Zpower_and_PSVs_v4.m
params.use_gamma_mod_chans = [true];%, false]; % [true, false] means all channels
% the following are not implemented & doing all instead
% target_encodingIDs = 1;
% target_correctness = 1;
% images_to_process = {}; %P046

% for ES
params.freq_min = min(params.ES_freq_band);
params.freq_max = max(params.ES_freq_band);

% between trials similarity
params.btwn_trial_type = 'EMS'; % default
params.comp_options = {'corr BT BI', 'corr BT WI'}; % between trial comparison types (between image or within image)
% note that compute_ES_btwn_trials_v6.m is only equipped to handle EMS, EES.

%% write this file,params to the output_folder

% saves default to empty if not passed.
% if save_script% ~isempty(save_script)
%     write_current_script_to_destination(fullfile(params.output_folder), "get_parameters.m");
%     if ~isempty(custom_params)
%         save(fullfile(params.output_folder,'custom_params.mat'), custom_params)
%     end
% end

% probably don't need this if we are saving custom_params and this script,
% but could use this for loading after first script?
% if save_values && ~isempty(save_values)
%     save(fullfile(params.output_folder,'parameters.mat'), params)
% end

%% override again in case any of the variables defined after the first override need adjustment too.

if nargin >= 1 && ~isempty(custom_params)
    fnames = fieldnames(custom_params);
    for i = 1:numel(fnames)
        % If the custom field exists, override the default value.
        % (Extra fields not already in params produce a warning.)
        if isfield(params, fnames{i})
            params.(fnames{i}) = custom_params.(fnames{i});
        else
            warning('Unknown parameter: %s. Ignored.', fnames{i});
        end
    end
end

%% init parallel pool of workers

if isempty(gcp('nocreate')) && params.init_par_pool% Open a local parallel pool if none exists
    if params.hellbender
        parpool('local', params.num_wrkrs_hb)
    else
        parpool('local', params.num_wrkrs_local);
    end
end

end