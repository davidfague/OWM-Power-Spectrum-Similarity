%%PSS_post_processing

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
    compute_all_within_trial_similarities(PS_file)
    % 
    % % gather into 1 file per patient since calculated separately first go
    % % around.
    combined_save_file = combine_enc_matrices(PS_file);
    % combined_save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', 'all3Enc_ES.mat'));
    % 
    ES_file = matfile(combined_save_file);
    [EMS, EFS, ~] = compute_EMS_EFS_EES(ES_file, PS_file);

    label_table = PS_file.label_table;
    rows_with_nan = any(isnan(label_table.EMS_means), 2);
    label_table = label_table(~rows_with_nan,:);
    label_table = subset_table_by_enc_performance(label_table);
    EMS_mean_byenc = mean(label_table.EMS_means(:,:),1);
    EFS_mean_byenc = mean(label_table.EFS_means(:,:),1);

    % sort by anatomy
    % anatomies = string(label_table.anatomical_label(:));
    % unique_labels = unique(anatomies);
    % EMS_mean_byenc_per_label = cell(length(unique_labels), 2);  % Create a cell array to store each label and corresponding mean
    % EFS_mean_byenc_per_label = cell(length(unique_labels), 2);  % Create a cell array to store each label and corresponding mea
    % for i = 1:length(unique_labels)
    % current_label = unique_labels{i};  % Get the current anatomical label
    % 
    % % Find the rows that match the current anatomical label
    % subset_indices = anatomies == current_label;
    % 
    % % Compute the mean for EMS_means for the current label
    % EMS_mean_byenc = mean(label_table.EMS_means(subset_indices, :), 1);
    % EFS_mean_byenc = mean(label_table.EFS_means(subset_indices, :), 1);
    % 
    % % Store the label and corresponding mean in the results container
    % EMS_mean_byenc_per_label{i, 1} = current_label;
    % EMS_mean_byenc_per_label{i, 2} = EMS_mean_byenc;
    % 
    % EFS_mean_byenc_per_label{i, 1} = current_label;
    % EFS_mean_byenc_per_label{i, 2} = EFS_mean_byenc;
    % end

    % Compute EMS_mean_byenc for each unique anatomical label using groupsummary
    % label_table.anatomical_label = string(label_table.anatomical_label);
    label_table.anatomical_label = strcat(string(label_table.anatomical_label), "_", num2str(patient_ID));
    % EMS_mean_byenc_per_label = groupsummary(label_table, 'anatomical_label', 'mean', 'EMS_means');
    % sorted_table = sortrows(EMS_mean_byenc_per_label, 'mean_EMS_means', 'descend');
    % Compute EMS_mean_byenc and EFS_mean_byenc for each unique anatomical label using groupsummary
    EMS_EFS_mean_byenc_per_label = groupsummary(label_table, 'anatomical_label', 'mean', {'EMS_means', 'EFS_means'});
    EMS_EFS_mean_byenc_per_label.Average_EMS = mean(EMS_EFS_mean_byenc_per_label.mean_EMS_means, 2);
    EMS_EFS_mean_byenc_per_label.Average_EFS = mean(EMS_EFS_mean_byenc_per_label.mean_EFS_means, 2);
    EMS_EFS_mean_byenc_per_label.diff = EMS_EFS_mean_byenc_per_label.mean_EMS_means - EMS_EFS_mean_byenc_per_label.mean_EFS_means;
    EMS_EFS_mean_byenc_per_label.avg_diff = mean(EMS_EFS_mean_byenc_per_label.diff, 2);


    
    % Sort the table by the mean of EMS_means in descending order
    % sorted_table = sortrows(EMS_EFS_mean_byenc_per_label, 'Average_EMS', 'descend');
    sorted_table = sortrows(EMS_EFS_mean_byenc_per_label, 'avg_diff', 'descend');

        % Concatenate the current patient's table with the combined table
    combined_table = [combined_table; EMS_EFS_mean_byenc_per_label];

    % PS_file.all_all3Enc_wholeTrial_ES_matrix() % a slice is row_id, enc_window_ids, all_window_ids, encID)
end
% Sort the combined table by the average difference in descending order
sorted_combined_table = sortrows(combined_table, 'avg_diff', 'descend');

function [EMS, EFS, EES] = compute_EMS_EFS_EES(ES_file, PS_file)
    [fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs, non_selection_win_IDs, all_win_IDs] = get_window_IDs(); % won't change each function call
    EMS = get_enc_similarity(ES_file, maint_win_IDs);
    EFS = get_enc_similarity(ES_file, fixation_win_IDs);
    EMS_means = compute_mean_ES(EMS);
    EFS_means = compute_mean_ES(EFS);
    % size(EMS_means) % 12180           1           1           3
    % size(EFS_means) % 12180           1           1           3
    % size(PS_file.label_table) % 12180           7

    % Squeeze singleton dimensions
    % EMS_means = squeeze(EMS_means);  % This will turn it into a 12180x3 matrix
    % EFS_means = squeeze(EFS_means);  % Same here

    % PS_file = matfile(fullfile(PS_file.Properties.Source), 'Writable', true);  % Ensure the file is writable

        % Load label_table from PS_file
    label_table = PS_file.label_table;
    

    % Update label_table in memory
    label_table.EMS_means = squeeze(EMS_means);
    label_table.EFS_means = squeeze(EFS_means);
    
    % Save the updated label_table back to PS_file
    % PS_file.label_table = label_table;
    save(PS_file.Properties.Source, 'label_table', '-append')
    EES = [];
    
end

function [label_table] = subset_table_by_enc_performance(label_table)
imageIDs = 1:9;
encIDs = 1:3;
target_encID = 3;

% size(label_table.encoding_correctness)
% rows = label_table.enc_correctness(:,target_encID) == 1; % rows where 3rd item was correct
rows = sum(label_table.encoding_correctness(:,:),2) == 3;
indices = rows;

label_table = label_table(indices,:);


end

function [mean_EFS_for_region, mean_EMS_for_region] = analyze_mean_ES_by_anat(label_table)
values = categorical(label_table.anatomical_label);
    for i=1:length(values)
        % Value to search for
        val = values(i);
        
        % Logical array where the condition is true
        rows = label_table.anatomical_label == val;
        
        % Get the indices of rows where condition is true
        indices = find(rows);
        
        % calculate mean EFS each region, mean EMS each region
        size(table(indices.EFS_means))
        mean_EFS_for_region = mean(table(indices).EFS_means);
    end
end

function means = compute_mean_ES(ES_matrix)
    % ES_means = zeros(size(ES_matrix, 1), 3);
    % for i=1:size(ES_matrix, 1) % channel-trial combinations
    mean_over_time = mean(ES_matrix, 3);

    means = mean(mean_over_time, 2); % mean across select_window_IDs then mean across encoding

    % means should be channel-trial x 3encodings
end

function similarity_subset = get_enc_similarity(PS_file, select_window_ids)
    similarity_subset = zeros(size(PS_file.all3_ES_matrix, 1), size(PS_file.all3_ES_matrix, 2), length(select_window_ids), 3);
    for i = 1:size(PS_file.all3_ES_matrix, 1)
        similarity_subset(i, :,:,:) = PS_file.all3_ES_matrix(i, :, select_window_ids, :); % each encoding vs part of trial
    end
end



function [] = compute_all_within_trial_similarities(data_file)


    [~, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, ~, non_selection_win_IDs, ~] = get_window_IDs(); % won't change each function call
    
    % data_file.all_windowed_mean_PS_vectors: size = all_window_IDs, frequencies, channels x trials 
    % label_table = data_file.label_table;

    % whole_trial_ES_matrix = compute_similarity_matrix(mean_PS_vectors, enc1_win_IDs, all_win_IDs);
    fprintf('\n')

    length_table = data_file.label_table;
    length_table = length(length_table.patient_ID);
    mean_PS_vectors = data_file.all_windowed_mean_PS_vectors(non_selection_win_IDs,:,:);
    encID_order = [1,3,2];
    all3_ES_matrix = zeros(length_table, length(enc1_win_IDs), length(non_selection_win_IDs), 3);
    save_file = fullfile(strrep(data_file.Properties.Source, '.mat', sprintf('all3_ES.mat')));
    for encID_idx = 1:3
        encID = encID_order(encID_idx);  % This will access enc1, enc3, enc2
        % save_file = fullfile(strrep(data_file.Properties.Source, '.mat', sprintf('_ES%d.mat', encID)));
        if exist(save_file, 'file')
            fprintf('%s already exists. Skipping computation.\n', save_file);
            return;
        else
            fprintf('computing %s\n', save_file);
        end

        % Determine the appropriate window IDs based on encID
        % Assign the correct window IDs based on encID
        switch encID
            case 1
                enc_win_IDs = enc1_win_IDs;
            case 2
                enc_win_IDs = enc2_win_IDs;
            case 3
                enc_win_IDs = enc3_win_IDs;
        end

        % temp_matrix = zeros(length_table, length(enc_win_IDs), length(non_selection_win_IDs), 3);
        
        % Compute similarity for all patients for this encoding ID
        parfor j = 1:length_table
            all3_ES_matrix(j, :, :, encID) = compute_similarity_matrix(mean_PS_vectors(:,:,j), enc_win_IDs, non_selection_win_IDs);
        end

        % save(save_file, 'temp_matrix', '-v7.3')

    % If the matrix already exists, skip execution

    % if isprop(save_file, 'all_all3Enc_wholeTrial_ES_matrix')
    %     disp('all_all3Enc_wholeTrial_ES_matrix already exists. Skipping computation.');
    %     return;
    % else
    %     fprintf('computing %s all_all3Enc_wholeTrial_ES_matrix', string(data_file));
    % end
            %all_all3Enc_wholeTrial_ES_matrix = zeros(length(label_table.patient_ID) ,size(whole_trial_ES_matrix, 1),size(whole_trial_ES_matrix, 2), 3);
        %        parfor i=1:length(label_table.patient_ID) % every channel-trial combination    
        %     % Temporary variable to hold the similarity matrices for this iteration
        %     all3Enc_wholeTrial_ES_matrix = zeros(size(whole_trial_ES_matrix, 1), size(whole_trial_ES_matrix, 2), 3);
        %     % Compute similarity matrices for each encoding period
        %     all3Enc_wholeTrial_ES_matrix(:, :, 1) = compute_similarity_matrix(mean_PS_vectors, enc1_win_IDs, all_win_IDs);
        %     all3Enc_wholeTrial_ES_matrix(:, :, 2) = compute_similarity_matrix(mean_PS_vectors, enc2_win_IDs, all_win_IDs);
        %     all3Enc_wholeTrial_ES_matrix(:, :, 3) = compute_similarity_matrix(mean_PS_vectors, enc3_win_IDs, all_win_IDs);
        % 
        %     % all3Enc_EMS_matrix = all3Enc_wholeTrial_ES_matrix(:,maint_win_IDs,:);
        %     % all3Enc_EFS_matrix = all3Enc_wholeTrial_ES_matrix(:,fixation_win_IDs,:);
        % 
        %     % all3Enc_EES_matrix = all3Enc_wholeTrial_ES_matrix(:,fixation_win_IDs,:);
        % 
        %     % Assign the temporary result to the main matrix
        %     all_all3Enc_wholeTrial_ES_matrix(i, :, :, :) = all3Enc_wholeTrial_ES_matrix;
        % end
    end

    save(save_file, 'all3_ES_matrix', '-v7.3')

end

function [combined_save_file] = combine_enc_matrices(data_file)

    % Precompute window IDs (won't change for each function call)
    % [fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs, non_selection_win_IDs, all_win_IDs] = get_window_IDs();
    
    % Output save file for the combined matrix
    combined_save_file = fullfile(strrep(data_file.Properties.Source, '.mat', 'all3Enc_ES.mat'));
    
    if exist(combined_save_file, 'file')
        fprintf('%s already exists. Skipping computation.\n', combined_save_file);
        return;
    end
    
    % Prepare a matrix to hold the combined similarities for each encoding period
    all3_ES_matrix = [];
    
    for encID = 1:3
        save_file = fullfile(strrep(data_file.Properties.Source, '.mat', sprintf('_ES%d.mat', encID)));
        
        if exist(save_file, 'file')
            fprintf('Loading %s\n', save_file);
            load(save_file, 'temp_matrix');
        else
            fprintf('ES%d file not found. Skipping this encoding.\n', encID);
            continue;
        end
        
        % Add a new dimension for each encoding matrix
        if isempty(all3_ES_matrix)
            % Initialize the combined matrix with the first encoding
            all3_ES_matrix = zeros(size(temp_matrix, 1), size(temp_matrix, 2), size(temp_matrix, 3), 3);
        end
        
        % Store the temp_matrix in the appropriate encoding ID slice
        all3_ES_matrix(:,:,:,encID) = temp_matrix;
    end
    
    % Save the combined matrix into a new file
    fprintf('Saving combined encoding matrices to %s\n', combined_save_file);
    save(combined_save_file, 'all3_ES_matrix', '-v7.3');
    
end
