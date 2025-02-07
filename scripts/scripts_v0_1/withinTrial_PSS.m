cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
addpath 'Raw Data Storage'               
addpath 'subfunctions'
%% TODO include other encodings.
% Define the folder and the criteria
folder_to_analyze = 'D:\Power Spectrum Similarity\output\201906 gammamod centralgyrus anyitem enc3 correct';
post_processing_folder = strrep(folder_to_analyze, 'output','output1');
if ~exist(post_processing_folder, 'dir')
    mkdir(post_processing_folder)
end
    %'201906 gammamod frontal anyitem enc3 correct';
%'D:\Power Spectrum Similarity\output\201906 gammamod centralgyrus anyitem enc3 correct';
%'C:\Users\david\Downloads\Power Spectrum Similarity\output\201906 gammamod centralgyrus anyitem enc3 correct';%201907 gammamod Amyg P046 enc3 correct';

patients_to_analyze = [201907];%[201906];
images_to_analyze = {};%{'P046'};
brain_locations_to_analyze = {};%{'central gyrus'};%{'amyg'};
channels_to_analyze = {};%{'76'};
trials_to_analyze = {};  % Empty means no filtering on this
encID_to_analyze = 1;
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
window_IDs2 = true(size(time));
window_IDs2 = find_valid_window_IDs_from_ntimes_logical_array(window_IDs2, window_start_times, window_end_times);
% Define the time windows for each encoding period
stim1_start = 1000; stim1_end = 1500; % stim 1 time window
stim2_start = 1500; stim2_end = 2000; % stim 2 time window
stim3_start = 2000; stim3_end = 2500; % stim 3 time window

%% new windows to compare against
fixation_times = (time >= 0 & time < 1000);
fixation_window_IDs = find_valid_window_IDs_from_ntimes_logical_array(fixation_times, window_start_times, window_end_times);

maintenance_times = (time >= 2500 & time < 6500);
maintenance_window_IDs = find_valid_window_IDs_from_ntimes_logical_array(maintenance_times, window_start_times, window_end_times);

%% Get a list of all files and folders in the specified directory
fileList = dir(folder_to_analyze);

% Filter out non-directories and exclude '.' and '..'
dirFlags = [fileList.isdir];  % Get a logical array for directories
subFolders = fileList(dirFlags);  % Keep only the directories
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % Exclude '.' and '..'

% Filter the list of folders
filtered_folders = filter_subfolders(subFolders, patients_to_analyze, brain_locations_to_analyze, channels_to_analyze, gamma_chans_to_analyze, images_to_analyze, trials_to_analyze, encID_to_analyze, encCorrectness_to_analyze);

% Display or use the filtered folder names
disp('Filtered subfolders that meet the criteria:')
for k = 1:length(filtered_folders)
    disp(filtered_folders(k).name);
end

if length(filtered_folders) < 1
    error('no subfolders to analyze')
end

%% do computations on certain parts of the mean_PSD_vectors or ES_matrix.
%% get fixation mean-power-spectrum vectors and compute similarity with encoding
averaging_data_folder = fullfile(post_processing_folder, 'averages across trials');
if ~exist(averaging_data_folder, 'dir')
    mkdir(averaging_data_folder)
end

averaged_data_savefile = fullfile(averaging_data_folder, 'across_trials.mat');
mean_fixation_similarity = cell(length(filtered_folders));
mean_maintenance_similarity = cell(length(filtered_folders));
MF_diff_mean_similarity = zeros(1, length(filtered_folders));
data_file_name = 'wholeTrial_ES.mat'; %'wholeTrial_ES_vectors.mat'
if ~exist(averaged_data_savefile, 'file')
    for file_id = 1:length(filtered_folders)
        % dir(fullfile(folder_to_analyze, filtered_folders(file_id).name));
        data_file_location = fullfile(folder_to_analyze, filtered_folders(file_id).name);
        data_file = fullfile(data_file_location, data_file_name);
        data_file = load(data_file);
        % data_file.similarity_matrix
        % data_file.window1_mean_PS_vectors
        % data_file.window2_mean_PS_vectors
    
        % make sure the window IDs are correct (only needed if windowing
        % changed or not using whole trial
        % [~, indices] = ismember(fixation_window_IDs, data_file.window_IDs1);  % use window_IDs2 bc it should be the whole trial
        % indices(indices == 0) = [];  % Remove zeros for items not found
        % fixation_window_ids_to_use = data_file.window_IDs2(1, indices);
        % if length(fixation_window_ids_to_use) ~= length(fixation_window_IDs)
        %     error('the desired windows are not all among the data.') % shouldn't happen if window_IDs2 is whole trial and using same windowing as before
        % end
    
        % get encoding mean power spectrum vectors for computing similarities
        % (no longer needed since we can subset existing ES matrix)
        % encoding_mean_PS_vectors = data_file.window1_mean_PS_vectors;
        %% compute encoding-fixation similarity
        % fixation_mean_PS_vectors = data_file.window2_mean_PS_vectors(fixation_window_IDs, :);
        % EF_similarity_matrix = compute_similarity_matrix_from_power_vecs(encoding_mean_PS_vectors, fixation_mean_PS_vectors);
        EFS_matrix = data_file.similarity_matrix(:,fixation_window_IDs); % alternative that subsets ES_matrix
        EFS_file = fullfile(data_file_location, 'EFS.mat');
        save(EFS_file, 'EFS_matrix')%,'-append')
    
        %% compute encoding-maintenance similarity
        % maintenance_mean_PS_vectors = data_file.window2_mean_PS_vectors(maintenance_window_IDs, :);
        % EM_similarity_matrix = compute_similarity_matrix_from_power_vecs(encoding_mean_PS_vectors, maintenance_mean_PS_vectors);
        EMS_matrix = data_file.similarity_matrix(:,maintenance_window_IDs); % alternative that subsets ES_matrix
        EMS_file = fullfile(data_file_location, 'EMS.mat');
        save(EMS_file, 'EMS_matrix')%,'-append')
        % alternative that subsets ES_matrix
    
        %% average similarity in region and PSD/
        % squash encoding time dimension
        average_EFS_vector = mean(EFS_matrix, 1);
        average_EMS_vector = mean(EMS_matrix, 1);
        % squash other time dimension as well
        average_EFS = mean(average_EFS_vector);
        average_EMS = mean(average_EMS_vector);
        % save extra computations
        save(EMS_file, "average_EMS_vector", "average_EMS", "-append")
        save(EFS_file, "average_EFS_vector", "average_EFS", "-append")
    
        figure_name = fullfile(data_file_location, 'avg_similarity.fig');
        if ~exist(figure_name, "file")
            EFS_file = matfile(EFS_file);
            EMS_file = matfile(EMS_file);
            fig = generate_avg_EMS_vs_EFS_figure(EFS_file, EMS_file, data_file_location);
            savefig(fig, figure_name);
            close(fig);
        end
    
        %% gather these averages; save these plots.
        mean_maintenance_similarity{file_id} = average_EMS;
        mean_fixation_similarity{file_id} = average_EFS;
        MF_diff_mean_similarity(file_id) = average_EMS - average_EFS;
    
    
        %% gather average withinTrial data across trials
        % gather shared non_trial_names; trial_folders for each non_trial_name
        % pool data from trial_folders for this non_trial_name
        % compute average ES for this channel-item-etc.
        if file_id == 1
            all_ES_matrices = zeros(size(data_file.similarity_matrix,1),size(data_file.similarity_matrix,2), length(filtered_folders));
            %% check size of data_file.window1_mean_PS_vectors and size of all_window1_mean_PS_vectors
            all_window1_mean_PS_vectors = zeros(size(data_file.window1_mean_PS_vectors, 1), size(data_file.window1_mean_PS_vectors, 2), length(filtered_folders)); % window 1 is encoding
            all_window2_mean_PS_vectors = zeros(size(data_file.window2_mean_PS_vectors,1), size(data_file.window2_mean_PS_vectors, 2), length(filtered_folders)); % window 2 is whole trial
        end
        all_ES_matrices(:,:, file_id) = data_file.similarity_matrix;
        all_window1_mean_PS_vectors(:,:, file_id) = data_file.window1_mean_PS_vectors;
        all_window2_mean_PS_vectors(:,:, file_id) = data_file.window2_mean_PS_vectors;
    end
    avg_ES = mean(all_ES_matrices, 3);
    avg_window1_mean_PS_vectors = mean(all_window1_mean_PS_vectors, 3); % check dimension hhere
    avg_window2_mean_PS_vectors = mean(all_window2_mean_PS_vectors, 3);
    save(averaged_data_savefile, "avg_ES", "avg_window1_mean_PS_vectors", "avg_window2_mean_PS_vectors")
    clear all_ES_matrices all_window1_mean_PS_vectors all_window2_mean_PS_vectors
else
    load(averaged_data_savefile)
end

%% save distribution of avg EMS- EFS across these trials
EMS_EFS_avgs_folder = fullfile(post_processing_folder, 'EMS_vs_EFS');
if ~exist(EMS_EFS_avgs_folder, 'dir')
    mkdir(EMS_EFS_avgs_folder)
end
title_from_folder = strsplit(folder_to_analyze, '\');
title_from_folder = title_from_folder{end};
EMS_EFS_avgs_savefile = fullfile(EMS_EFS_avgs_folder, 'EMS_vs_EFS.mat');
save(EMS_EFS_avgs_savefile, 'mean_maintenance_similarity', 'mean_fixation_similarity', 'MF_diff_mean_similarity')
clear EMS_EFS_avgs_savefile mean_maintenance_similarity mean_fixation_similarity
%% plot distribution of MF_diff_mean_similarity
fig = figure('WindowState','maximized');%, 'Visible', 'off');
histogram(MF_diff_mean_similarity)
title(['Maint-specific similarity, Distribution of mean(EMS)-mean(EFS) across ' title_from_folder]);
xlabel('Value');
ylabel('Frequency');
figure_file = fullfile(EMS_EFS_avgs_folder, 'distribution_of_EMSavg-EFSavg_each_trial.fig');
savefig(fig, figure_file)
% close(fig)
clear EMS_EFS_avgs_folder MF_diff_mean_similarity
%% save and plot average whole-within-Trial ES across these trials
avg_ES_folder = fullfile(post_processing_folder, 'avg_wholeTrial_ES');
if ~exist(avg_ES_folder, 'dir')
    mkdir(avg_ES_folder)
end
avg_ES_savefile = fullfile(avg_ES_folder, 'avg_wholeTrial_ES.mat');
save(avg_ES_savefile, 'avg_ES', 'avg_window1_mean_PS_vectors', 'avg_window2_mean_PS_vectors')
fig = generate_avg_ES_figure(avg_ES, title_from_folder);
fig_savename = fullfile(avg_ES_folder, 'avg_wholeTrial_ES.fig');
savefig(fig, fig_savename)
% close(fig)

% save and plot average mean PSD vectors across these trials
fig = generate_ES_vecs_figure(matfile(avg_ES_savefile), post_processing_folder);
fig_savename = fullfile(avg_ES_folder, 'avg_wholeTrial_ES_vectors.fig');
savefig(fig, fig_savename)
% close(fig)
clear 'avg_ES' 'avg_window1_mean_PS_vectors' 'avg_window2_mean_PS_vectors' 'fig_savename' 'avg_ES_folder' 'avg_ES_savefile'
%%
function fig = generate_avg_ES_figure(avg_ES, folder_to_analyze)
    % can turn this into function and clean up
    fig = figure('WindowState','maximized');%, 'Visible', 'off');%, 'Visible', 'off');
    title_str = strrep(folder_to_analyze, '_', ' ');
    title_str = strcat('Average Within-Trial Whole-Trial Encoding Similarity', title_str);
    title(title_str, 'Interpreter', 'none')
    imagesc(avg_ES)
    c = colorbar;
    c.Label.String = 'Similarity (rho)';
    clim([0 1])
    title(title_str);
    xlabel('Whole Trial Window ID');
    ylabel('Encoding Window ID');
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

function fig = generate_ES_vecs_figure(ES_file, folder_name)
    fig = figure('WindowState','maximized');%, 'Visible', 'off');%,'Visible', 'off');
    
    title_str = strrep(folder_name, '_', ' ');
    title_str = strcat("Mean Z Power By-Frequency Vectors ", title_str);

    % Calculate the relative widths based on the data sizes
    encoding_data_size = size(ES_file.avg_window1_mean_PS_vectors, 1); % Encoding x-axis size (40)
    trial_data_size = size(ES_file.avg_window2_mean_PS_vectors, 1); % Whole trial x-axis size (891)
    total_size = encoding_data_size + trial_data_size;
 
    encoding_width = encoding_data_size / total_size; % Proportional width for encoding subplot
    trial_width = trial_data_size / total_size; % Proportional width for whole trial subplot
    
    % First subplot: Display window1_mean_PS_vectors (Item Encoding)
    subplot(2, 1, 1); % Create a 2x1 grid, first plot
    imagesc(ES_file.avg_window1_mean_PS_vectors');
    c=colorbar;
    c.Label.String = 'Z-scored Power';
    title('Item Encoding', 'Interpreter', 'none');
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
    imagesc(ES_file.avg_window2_mean_PS_vectors');
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
    % Get the upper limit of the y-axis for label placement
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
    sgtitle(title_str, 'Interpreter', 'none'); % Set super title for the figure
    
    % Optionally, set colormap for both subplots
    % colormap('hot'); % Apply a colormap
end

function fig = generate_avg_EMS_vs_EFS_figure(EFS_file, EMS_file, folder_name)
    % can turn this into function and clean up
    fig = figure('WindowState','maximized', 'Visible', 'off');%, 'Visible', 'off');
    title_str = strrep(folder_name, '_', ' ');
    title_str = strcat('Average Within-Trial Encoding Similarity ', title_str);
    title(title_str)
    plot(EFS_file.average_EFS_vector, 'color','r')
    hold on
    plot(EMS_file.average_EMS_vector, 'color','b')
    yline(EFS_file.average_EFS, 'color','red')
    yline(EMS_file.average_EMS, 'color','blue')
    ylabel('Similarity (rho)')
    xlabel('time window index')
    legend('enc-fixation', 'enc-maint', 'average EFS','average EMS')
end

%% Now create example of mix-matching trial similarity (i.e. random pairing of trial1 enc with trial2 maint)

%% also create plotting for new similarity matrices.

%% generating clusters for test
%% null distribution is

%% motor regionsf