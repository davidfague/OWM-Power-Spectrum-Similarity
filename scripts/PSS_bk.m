% todo check power computation with another method; fix zbaseline for single trial
%% consider renaming window1_mean_PS_vectors; reorganizing saved folders
%% Power Spectrum Similarity
%cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
addpath 'Raw Data Storage'               
addpath 'subfunctions'
%% Specify data & subset - regions, gamma/non-gamma channels, trials, performance, items
patient_IDs = [201907 ];
use_gamma_mod_chans = [true];%, false];
images_to_process = {}; %P046
brain_anatomies_to_process = {'frontal'};
target_encodingIDs = 3;
target_correctness = 1;
output_folder = 'D:\Power Spectrum Similarity\output\allpatients gammamod allregions allitem enc3 correct';
%'C:\Users\david\Downloads\Power Spectrum Similarity\output\201907 gammamod hipp anyitem enc3 correct';
%% Similarity processing
ES_prefix = 'wholeTrial';
time = 1:9001;
% 100ms time window with step of 10ms
wt = length(time);
num_windows = floor((wt - 100) / 10) + 1; % Correct calculation for number of windows
% Precompute window start and end times
window_start_times = time(1:10:(num_windows - 1) * 10 + 1); % Start times
clear num_windows wt
window_end_times = time(window_start_times + 99);            % End times

% access window start and end: (:, window_id)
% windows for 2nd iteration (whole trial windows)
window_IDs2 = true(size(time));
window_IDs2 = find_valid_window_IDs_from_ntimes_logical_array(window_IDs2, window_start_times, window_end_times);
% Define the time windows for each encoding period
stim1_start = 1000; stim1_end = 1500; % stim 1 time window
stim2_start = 1500; stim2_end = 2000; % stim 2 time window
stim3_start = 2000; stim3_end = 2500; % stim 3 time window

ES_freq_band = 1:40;

save_mean_PS_vectors = true;

%% figure constants
colorbar_limits = [0, 1];

%% prepare to load data using matfile for memory efficiency

% Preallocate cell arrays for data_matrix and data_info matfile objects
patient_data_files = cell(1, length(patient_IDs));          % store data files
patient_info_files = cell(1, length(patient_IDs));          % store info files
patient_images = cell(1, length(patient_IDs));              % 9 images per patient
patient_enc_to_imageID = cell(1, length(patient_IDs));      % To store encoding itemID      for each each trial, enc, patient
patient_enc_correctness = cell(1, length(patient_IDs));     % To store encoding correctness for each trial, enc, patient
patient_is_gamma_channels = cell(1, length(patient_IDs));   % To store gamma modulation     for each channel, patient
patient_original_channel_IDs = cell(1, length(patient_IDs));% store original channel IDs per patient
patient_original_trial_IDs  = cell(1, length(patient_IDs)); % store original trial IDs per patient
% Iterate over patient IDs and store matfile objects in the cell arrays
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);  % Assuming patient_IDs is a cell array of strings
    % Store matfile objects in the corresponding cell array
    patient_data_files{idx} = matfile(fullfile('Raw Data Storage', ['D_OWM_t_bipolar_', num2str(patient_ID), '.mat']));
    patient_info_files{idx} = matfile(fullfile('Raw Data Storage', ['OWM_trialinfo_', num2str(patient_ID), '.mat']));

    % load patient images
    patient_images{idx} = patient_info_files{idx}.C; % Load the full variable

    % Load Cond_item (nTrial x 3 enc) and Cond_performance (nTrial x 3 enc) for each patient
    patient_enc_to_imageID{idx} = patient_info_files{idx}.Cond_item;
    patient_enc_correctness{idx} = patient_info_files{idx}.Cond_performance;

    % load channel gamma modulation
    gamma_data = load(fullfile('Raw Data Storage', ['gammachans2sd_alltrials_', num2str(patient_ID), '.mat']));
    num_channels = size(patient_data_files{1}.labelsanatbkedit, 1);         
    patient_is_gamma_channels{idx} = bool_mask_array(gamma_data.sigchans2, num_channels);

    % store original IDs to subset and carry later
    patient_original_channel_IDs{idx} = 1:num_channels;
    num_trials = size(patient_enc_correctness{idx}, 1);
    patient_original_trial_IDs{idx} = 1:num_trials;
    
end
clear idx patient_ID gamma_data num_channels num_trials

%% subset data by gamma mod & anat channels;
patient_channels_to_process = cell(length(patient_IDs));
patient_channel_brain_locations_to_process = cell(length(patient_IDs));
for idx = 1:length(patient_IDs)
    % choose channels - subset original channel IDs
    % choose channels by gamma modulation - subset original channel IDs
    patient_channels_to_process{idx} = subset_channels_by_gamma_modulation(patient_original_channel_IDs{idx}, use_gamma_mod_chans, patient_is_gamma_channels{idx});

    % choose channels by brain location %- label by location.
    patient_channel_brain_locations = patient_data_files{idx}.labelsanatbkedit; % matfile objects does not support partial loading of tables
    patient_channel_brain_locations = patient_channel_brain_locations.anatmacro1; % matfile objects does not support partial loading of tables
    patient_channels_to_process{idx} = subset_channels_by_brain_location(patient_channels_to_process{idx}, brain_anatomies_to_process, patient_channel_brain_locations);
    patient_channel_brain_locations_to_process{idx} = patient_channel_brain_locations(patient_channels_to_process{idx});
end
clear idx patient_original_channel_IDs use_gamma_mod_chans

%% subset trialIDs by image, encCorrectness and encID
patient_image_ids_to_keep = cell(length(patient_IDs));
patient_image_to_trialID_encodingID_encodingCorrectness = cell(length(patient_IDs));
for idx = 1:length(patient_IDs)
    % choose images - subset trials by images; prepare image labels
    [patient_image_ids_to_keep{idx}, patient_image_to_trialID_encodingID_encodingCorrectness{idx}] = subset_trials_by_image_occurence( ...
        patient_original_trial_IDs{idx}, images_to_process, patient_enc_to_imageID{idx}, patient_images{idx}, patient_enc_correctness{idx});

    % subset by target encoding correctness and encodingID
    patient_image_to_trialID_encodingID_encodingCorrectness{idx} = subset_by_encodingCorrectness_and_ID(...
        patient_image_ids_to_keep{idx}, patient_image_to_trialID_encodingID_encodingCorrectness{idx}, target_correctness, target_encodingIDs);
end

clear idx images_to_process brain_anatomies_to_process target_encodingIDs target_correctness
%% compute power for data subset

% data_subset = patient_data_files{patient_idx}.D_OWM_t(patient_channels_to_process{idx}, :,patient_image_to_trialID_encodingID_encodingCorrectness{patient_idx}{:,1});

folder_list = {};

% ID corresponds to index in original matrix; idx corresponds to index in subset
for patient_idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(patient_idx);

    % preallocate new matrix for this subset
    num_channels = length(patient_channels_to_process{patient_idx});
    % num_trials = length(patient_image_to_trialID_encodingID_encodingCorrectness{patient_idx}{:,1});
    num_timepoints = 9001; % 9001 time points

    for channel_idx = 1:length(patient_channels_to_process{patient_idx})
        channel_ID = patient_channels_to_process{patient_idx}(channel_idx);
        brain_location = patient_channel_brain_locations_to_process{patient_idx}{channel_idx};
        brain_location = strrep(brain_location, ' ','');
        channel_is_gamma = patient_is_gamma_channels{patient_idx}(channel_ID);
        for image_idx = 1:length(patient_image_ids_to_keep{patient_idx})
            image_name = patient_images{patient_idx}{patient_image_ids_to_keep{patient_idx}};
            image_name = strrep(image_name, '.', '');
            for trial_idx = 1:length(patient_image_to_trialID_encodingID_encodingCorrectness{patient_idx}{image_idx,1})
                trial_ID = patient_image_to_trialID_encodingID_encodingCorrectness{patient_idx}{image_idx,1}(trial_idx);
                enc_ID = patient_image_to_trialID_encodingID_encodingCorrectness{patient_idx}{image_idx,2}(trial_idx);
                enc_Correctness = patient_image_to_trialID_encodingID_encodingCorrectness{patient_idx}{image_idx,3}(trial_idx);

                % Load the specific data slice for this channel and trial
                data_subset = patient_data_files{patient_idx}.D_OWM_t(channel_ID, :, trial_ID);

                % check for nans - skip this data there are nans
                if any(isnan(data_subset(:)))
                    % Skip the current iteration if NaNs are found
                    continue;
                end
                                % Format the folder name
                folder_name = strcat('Patient', num2str(patient_ID), ...
                                     '_Loc', brain_location, ...
                                     '_Chan', num2str(channel_ID), ...
                                     '_ChanisGam', num2str(channel_is_gamma), ...
                                     '_Image', image_name, ...
                                     '_Trial', num2str(trial_ID), ...
                                     '_EncId', num2str(enc_ID), ...
                                     '_Correct', num2str(enc_Correctness));

                % Create the full path for the new folder
                data_subset_folder = fullfile(output_folder, folder_name);

                folder_list{end+1} = data_subset_folder;
                
                % Check if the folder already exists, if not, create it
                if ~exist(data_subset_folder, 'dir')
                    mkdir(data_subset_folder);
                end

                info_file = fullfile(data_subset_folder, 'info.mat');
                if ~exist(info_file, "file")
                    save(info_file, "patient_ID", "brain_location", "channel_ID", "channel_is_gamma", "image_name", "trial_ID", "enc_ID", "enc_Correctness")
                end

                data_subset_file = fullfile(data_subset_folder, 'raw_data.mat');
                if ~exist(data_subset_file, "file")
                    save(data_subset_file, "data_subset")
                end

                power_file = fullfile(data_subset_folder, 'Zpower.mat');
                if ~exist(power_file, "file")
                    Zpower = compute_Zpower(data_subset);
                    save(power_file, "Zpower")
                end

                ES_file = fullfile(data_subset_folder, strcat(ES_prefix, '_ES.mat'));
                if ~exist(ES_file, "file")
                    power_file_data = load(power_file);
                    % get ntimes logical array corresponding to this
                    % encoding timeframe
                    window_IDs1 = get_encoding_times_from_enc_ID(enc_ID, time, stim1_start, stim1_end, stim2_start, stim2_end, stim3_start, stim3_end);
                    % get window IDs from logical array
                    window_IDs1 = find_valid_window_IDs_from_ntimes_logical_array(window_IDs1, window_start_times, window_end_times);
                    % now have 2 logical arrays of size nTimes:
                    % window_IDs1, window_IDs2
                    %use logical arrays to subset window_IDs.

                    [similarity_matrix, window1_mean_PS_vectors, window2_mean_PS_vectors] = compute_ES(power_file_data.Zpower, window_start_times, window_end_times, window_IDs1, window_IDs2, ES_freq_band, save_mean_PS_vectors);
                    save(ES_file, "similarity_matrix", "window1_mean_PS_vectors", "window2_mean_PS_vectors", "window_IDs1", "window_IDs2", "ES_freq_band", "save_mean_PS_vectors")
               % elseif
                end

                % ES figure
                figure_name = fullfile(data_subset_folder, strcat(ES_prefix, '_ES.fig'));
                % if ~exist(figure_name, "file")
                    ES_file = matfile(ES_file);
                    fig = generate_ES_figure(ES_file, folder_name, colorbar_limits);
                    savefig(fig, figure_name);
                    close(fig)
                % end
                    
                % mean power vectors figure
                figure_name = fullfile(data_subset_folder, strcat(ES_prefix, '_ES_vectors.fig'));
                % if ~exist(figure_name, "file")
                    fig = generate_ES_vecs_figure(ES_file, folder_name, colorbar_limits);
                    savefig(fig, figure_name);
                    close(fig);
                % end
                % full_data_subset(channel_idx, :, trial_idx) = data_subset;
            end
        end
    end
end
clear data_subset

function fig = generate_ES_figure(ES_file, folder_name, colorbar_limits)
    fig = figure('WindowState','maximized', 'Visible', 'off');
    title_str = strrep(folder_name, '_', ' ');
    title_str = strcat('Within-Trial Encoding Similarity ', title_str);
    imagesc(ES_file.similarity_matrix);
    c = colorbar;
    c.Label.String = 'Similarity (rho)';
    clim(colorbar_limits)
    title(title_str);
    xlabel('Trial Time');
    ylabel('Encoding Time');
    % colormap('hot');

    % Plot the vertical lines
    index_adjust = -5; % subtracting 5 actually picks the window id that is centered at the desired time instead of the window id that begins at the desired time.
    xline(25+index_adjust, 'g', 'LineWidth', 2);    % baseline start
    xline(75+index_adjust, 'g', 'LineWidth', 2);    % baseline end
    xline(100+index_adjust, 'b', 'LineWidth', 2);   % fixation end
    xline(150+index_adjust, 'b', 'LineWidth', 2);   % enc1 end
    xline(200+index_adjust, 'b', 'LineWidth', 2);   % enc2 end
    xline(250+index_adjust, 'b', 'LineWidth', 2);   % enc3 end
    xline(650+index_adjust, 'b', 'LineWidth', 2);   % maintenance end
    label_loc = ylim;
    label_loc = label_loc(2);
    label_loc = label_loc + 5;
    % Add text labels at corresponding positions
    text(25+index_adjust, label_loc, 'baseline start', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'g');
    text(75+index_adjust, label_loc, 'baseline end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'g');
    text(100+index_adjust, label_loc, 'fixation end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(150+index_adjust, label_loc, 'enc1 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(200+index_adjust, label_loc, 'enc2 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(250+index_adjust, label_loc, 'enc3 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(650+index_adjust, label_loc, 'maint end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
end

function fig = generate_ES_vecs_figure(ES_file, folder_name, colorbar_limits)
    fig = figure('WindowState','maximized','Visible', 'off');
    
    title_str = strrep(folder_name, '_', ' ');
    title_str = strcat("Mean Z Power By-Frequency Vectors ", title_str);

    % Calculate the relative widths based on the data sizes
    encoding_data_size = size(ES_file.window1_mean_PS_vectors, 1); % Encoding x-axis size (40)
    trial_data_size = size(ES_file.window2_mean_PS_vectors, 1); % Whole trial x-axis size (891)
    total_size = encoding_data_size + trial_data_size;
 
    encoding_width = encoding_data_size / total_size; % Proportional width for encoding subplot
    trial_width = trial_data_size / total_size; % Proportional width for whole trial subplot
    
    % First subplot: Display window1_mean_PS_vectors (Item Encoding)
    subplot(2, 1, 1); % Create a 2x1 grid, first plot
    imagesc(ES_file.window1_mean_PS_vectors');
    c=colorbar;
    c.Label.String = 'Z-scored Power';
    title('Item Encoding');
    xlabel('Encoding Windows');
    ylabel('Frequency');
    
    % Capture color limits from first plot, increase range and use for
    % both plots
    limits_to_use = clim;
    limits_to_use(1) = limits_to_use(1) - 1;
    limits_to_use(2) = limits_to_use(2) + 1;
    clim(limits_to_use)

    % Adjust subplot width based on data size ratio
    ax1 = gca; % Get current axis
    pos1 = get(ax1, 'Position'); % Get current position
    pos1(3) = encoding_width * 0.85; % Adjust width (0.85 to keep margins)
    set(ax1, 'Position', pos1); % Apply new position

    % Second subplot: Display window2_mean_PS_vectors (Whole Trial)
    subplot(2, 1, 2); % Create a 2x1 grid, second plot
    imagesc(ES_file.window2_mean_PS_vectors');
    c=colorbar;
    clim(limits_to_use); % Apply same color limits as the first plot
    c.Label.String = 'Z-scored Power'; % Customize the label text as needed
    title('Whole Trial');
    xlabel('Trial Windows');
    ylabel('Frequency');
    
    % Plot the vertical lines
    index_adjust = -5; % subtracting 5 actually picks the window id that is centered at the desired time instead of the window id that begins at the desired time.
    xline(25+index_adjust, 'g', 'LineWidth', 2);    % baseline start
    xline(75+index_adjust, 'g', 'LineWidth', 2);    % baseline end
    xline(100+index_adjust, 'b', 'LineWidth', 2);   % fixation end
    xline(150+index_adjust, 'b', 'LineWidth', 2);   % enc1 end
    xline(200+index_adjust, 'b', 'LineWidth', 2);   % enc2 end
    xline(250+index_adjust, 'b', 'LineWidth', 2);   % enc3 end
    xline(650+index_adjust, 'b', 'LineWidth', 2);   % maintenance end
    label_loc = ylim;
    label_loc = label_loc(2);
    label_loc = label_loc + 5;
    % Add text labels at corresponding positions
    text(25+index_adjust, label_loc, 'baseline start', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'g');
    text(75+index_adjust, label_loc, 'baseline end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'g');
    text(100+index_adjust, label_loc, 'fixation end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(150+index_adjust, label_loc, 'enc1 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(200+index_adjust, label_loc, 'enc2 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(250+index_adjust, label_loc, 'enc3 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(650+index_adjust, label_loc, 'maint end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    % Adjust subplot width based on data size ratio
    ax2 = gca; % Get current axis
    pos2 = get(ax2, 'Position'); % Get current position
    pos2(3) = trial_width * 0.85; % Adjust width proportional to data size (0.85 to keep margins)
    set(ax2, 'Position', pos2); % Apply new position

    % Set overall figure title
    sgtitle(title_str); % Set super title for the figure
    
    % Optionally, set colormap for both subplots
    % colormap('hot'); % Apply a colormap
end

%% use folder_list yo do EM and EF analysis