%% within trial clustering of similarity
%% Power Spectrum Similarity
cd('D:\Power Spectrum Similarity')% cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
addpath 'Raw Data Storage'               
addpath 'subfunctions'


%% Specify large data subset to process - empty means all - regions, gamma/non-gamma channels, trials, performance, items
patient_IDs = [201907];%[201907 201908, 201903, 201905, 201906, 201901, 201910, 201915];
use_gamma_mod_chans = [true];%, false]; % [true, false] means all channels
brain_anatomies_to_process = {};
% target_encodingIDs = 1;
% target_correctness = 1;
% images_to_process = {}; %P046;

% output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc');
output_folder = 'D:\Power Spectrum Similarity\parallel output\allpatients gammamod allregions allitem enc1 correct';

% note remove pwd after moving data from pwd back to the correct folder.

%'parallel output\allpatients gammamod allregions allitem enc1 correct');
%'D:\Power Spectrum Similarity\parallel output\allpatients gammamod allregions allitem enc1 correct';
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

%% constants
[fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs] = get_window_IDs();

%% intialize parallel pool
% if isempty(gcp('nocreate'))
%     parpool('local', 64);%25); % Open a local parallel pool if none exists
% end
combined_table = table();
%% loop through patients
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    PS_file = matfile(strcat(fullfile(output_folder, num2str(patient_ID)), '.mat'));

    % % % compute correlations
    % compute_all_within_trial_similarities(PS_file)
    % 
    % % gather into 1 file per patient since calculated separately first go
    % % around.
    % combined_save_file = combine_enc_matrices(PS_file);
    combined_save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', 'all3Enc_ES.mat'));
    % 
    ES_file = matfile(combined_save_file);
    % [EMS, EFS, ~] = compute_EMS_EFS_EES(ES_file, PS_file);

    % label_table = PS_file.label_table;
    % rows_with_nan = any(isnan(label_table.EMS_means), 2);
    % label_table = label_table(~rows_with_nan,:);
    % % label_table = subset_table_by_enc_performance(label_table); % selects
    % % only all 3 correct
    % EMS_mean_byenc = mean(label_table.EMS_means(:,:),1);
    % EFS_mean_byenc = mean(label_table.EFS_means(:,:),1);
    % 
    % label_table.anatomical_label = strcat(string(label_table.anatomical_label), "_", num2str(patient_ID));

    compute_clusters(PS_file, ES_file)
end

function compute_clusters(PS_file, ES_file, patient_ID)
    % Load target label and file information
    target_anat = "R IFG (p Orbitalis)";
    
    % Step 1: Filter rows with NaNs in EMS_means
    label_table = PS_file.label_table;
    rows_with_nan = any(isnan(label_table.EMS_means), 2);
    label_table = label_table(~rows_with_nan, :);
    indices_without_nan = find(~rows_with_nan);
    
    % Step 2: Select rows with the target anatomical label
    label_table.anatomical_label = string(label_table.anatomical_label);
    rows_with_target_anat = label_table.anatomical_label == target_anat;
    indices_with_target_anat = find(rows_with_target_anat);
    
    if isempty(indices_with_target_anat)
        warning('No indices found with the target anatomical label.');
        return;
    end
    
    % Step 3: Initialize the results for selected indices
    selected_index = indices_with_target_anat(1); % Plot the first instance as an example
    
    % Step 4: Load data in chunks to optimize memory usage
    chunk_size = 1000; % Set a chunk size that fits into memory (adjust based on memory availability)
    num_chunks = ceil((max(indices_without_nan) - min(indices_without_nan) + 1) / chunk_size);
    
    % Initialize filtered matrices
    ES_matrix_filtered = [];
    PS_matrix_filtered = [];
    
    for chunk = 1:num_chunks
        % Define the chunk range
        chunk_start_idx = min(indices_without_nan) + (chunk - 1) * chunk_size;
        chunk_end_idx = min(chunk_start_idx + chunk_size - 1, max(indices_without_nan));
        
        % Step 5: Load the chunk of data from ES_file and PS_file
        full_ES_chunk = ES_file.all3_ES_matrix(chunk_start_idx:chunk_end_idx, :, :, :);
        full_PS_chunk = PS_file.all_windowed_mean_PS_vectors(:, :, chunk_start_idx:chunk_end_idx);
        
        % Step 6: Filter the chunk based on indices_without_nan
        valid_in_chunk = ismember(chunk_start_idx:chunk_end_idx, indices_without_nan);
        ES_matrix_filtered = cat(1, ES_matrix_filtered, full_ES_chunk(valid_in_chunk, :, :, :));
        PS_matrix_filtered = cat(3, PS_matrix_filtered, full_PS_chunk(:, :, valid_in_chunk));
    end
    
    % Step 7: Plot the results using the selected index
    main_fig = create_main_figure(ES_matrix_filtered, PS_matrix_filtered, selected_index, patient_ID);
    
    % Save or display the figure as needed
    savefig(main_fig, sprintf('Patient_%d_ClusterFigure.fig', patient_ID));
end

function main_fig = create_main_figure(ES_matrix, PS_matrix, index, patient_ID)
    % Create the main figure
    main_fig = figure('WindowState', 'maximized');
    
    % Create subplots for each component
    subplot(1, 2, 1); % First subplot for representational similarity
    plot_similarity(squeeze(ES_matrix(index, :, :, :)));
    
    subplot(1, 2, 2); % Second subplot for PS_matrix representation
    plot_representation(squeeze(PS_matrix(:, :, index)));
    
    % Add a title to the main figure
    sgtitle(sprintf('Patient %d - Cluster Visualization', patient_ID));
end

function plot_similarity(ES_matrix_slice)
    % Plot similarity heatmap with relevant annotations
    imagesc(ES_matrix_slice);
    colorbar;
    title('Encoding Similarity');
    xlabel('Encoding Windows');
    ylabel('Whole Trial Windows');
    hold on;

    % Plot vertical lines for various events
    add_event_lines();
end

function plot_representation(PS_matrix_slice)
    % Plot power spectrum heatmap with relevant annotations
    imagesc(PS_matrix_slice');
    colorbar;
    title('Representation Similarity (Power Spectrum)');
    xlabel('Whole Trial Window ID');
    ylabel('Frequency');
    hold on;

    % Plot vertical lines for various events
    add_event_lines();
end

% function compute_clusters(PS_file, ES_file)
% % size(PS_file.label_table) % (trials*channels, cols) (12180, 9)
% % size(ES_file.all3_ES_matrix) % (trials*channels, encodingWindows, wholeTrialWindows, encodings) (12180, 41, 640, 3)
% size(PS_file.all_windowed_mean_PS_vectors) % (wholeTrialWindows, frequencies, trials*channels) (891, 40, 12180)
% 
% target_anat = "R IFG (p Orbitalis)_201907";
% 
% label_table = PS_file.label_table;
% rows_with_nan = any(isnan(label_table.EMS_means), 2);
% label_table = label_table(~rows_with_nan,:);
% indices_without_nan = find(~rows_with_nan);
% ES_matrix = ES_file.all3_ES_matrix;
% ES_matrix = ES_matrix(indices_without_nan,:,:,:);
% 
% PS_matrix = PS_file.all_windowed_mean_PS_vectors(:,:,indices_without_nan);
% 
% % target anatomy label
% rows_with_target_anat = label_table.anatomical_label == target_anat;
% indices_with_target_anat = find(rows_with_target_anat);
% % label_table(label_table.anatomical_label == target_anat,:);
% 
% % plot
% plot_encoding_similarity(squeeze(ES_matrix(indices_with_target_anat(1),:,:,:)));
% 
% plot_representation(squeeze(PS_matrix(:,:,indices_with_target_anat(1))));
% 
% 
% % % label_table = subset_table_by_enc_performance(label_table); % selects
% % % only all 3 correct
% % EMS_mean_byenc = mean(label_table.EMS_means(:,:),1);
% % EFS_mean_byenc = mean(label_table.EFS_means(:,:),1);
% 
% label_table.anatomical_label = strcat(string(label_table.anatomical_label), "_", num2str(patient_ID));
% end

% function plot_similarity(ES_matrix)
% if length(size(ES_matrix)) > 3
%     warning("need to subset ES_matrix some more before plotting. ndim= %s", length(size(ES_matrix)))
% end
% 
% figure()
% imagesc(ES_matrix)
% end
% 
% function plot_representation(PS_matrix);
% if length(size(PS_matrix)) > 2
%     warning("need to subset PS_matrix some more before plotting. ndim= %s", length(size(ES_matrix)))
% end
% figure()
% imagesc(PS_matrix)
% end
% 
% function similarity_subset = get_enc_similarity(PS_file, select_window_ids)
%     similarity_subset = zeros(size(PS_file.all3_ES_matrix, 1), size(PS_file.all3_ES_matrix, 2), length(select_window_ids), 3);
%     for i = 1:size(PS_file.all3_ES_matrix, 1)
%         similarity_subset(i, :,:,:) = PS_file.all3_ES_matrix(i, :, select_window_ids, :); % each encoding vs part of trial
%     end
% end
