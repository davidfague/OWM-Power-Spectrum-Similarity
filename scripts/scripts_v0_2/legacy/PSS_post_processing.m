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

%% loop through patients - compute within_trial similarities
combined_table = table();
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

    label_table = PS_file.label_table;
    rows_with_nan = any(isnan(label_table.EMS_means), 2);
    label_table = label_table(~rows_with_nan,:);
    % label_table = subset_table_by_enc_performance(label_table); % selects
    % only all 3 correct
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
%% testing
    [EMS_EFS_mean_byenc_per_label] =  compute_stats_table(label_table);
%% testing
    % % Group by 'anatomical_label' and 'encoding_correctness'
    % EMS_EFS_mean_byenc_per_label = groupsummary(label_table, {'anatomical_label', 'encoding_correctness'}, 'mean', {'EMS_means', 'EFS_means'});
    % % Calculate average EMS and EFS across encoding correctness for each label
    % EMS_EFS_mean_byenc_per_label.Average_EMS = mean(EMS_EFS_mean_byenc_per_label.mean_EMS_means, 2);
    % EMS_EFS_mean_byenc_per_label.Average_EFS = mean(EMS_EFS_mean_byenc_per_label.mean_EFS_means, 2);
    % % Calculate the difference between EMS and EFS means
    % EMS_EFS_mean_byenc_per_label.diff = EMS_EFS_mean_byenc_per_label.mean_EMS_means - EMS_EFS_mean_byenc_per_label.mean_EFS_means;
    % % Calculate the average difference across encoding correctness for each label
    % EMS_EFS_mean_byenc_per_label.avg_diff = mean(EMS_EFS_mean_byenc_per_label.diff, 2);

    %% sort and combine tables
    % Sort the table by the mean of EMS_means in descending order
    % sorted_table = sortrows(EMS_EFS_mean_byenc_per_label, 'Average_EMS', 'descend');
    sorted_table = sortrows(EMS_EFS_mean_byenc_per_label, 'avg_diff', 'descend');
    % save("by_performance.mat", 'sorted_table')
    % sorted_table = sortrows(EMS_EFS_mean_byenc_per_label, 'anatomical_label', 'descend');
    [difference_table] = compute_diff_btwn_corr_incorr(sorted_table);
    difference_table = sortrows(difference_table, 'avg_diff_Difference', 'descend');
    % save("performance_summary.mat", 'difference_table')
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

function [EMS_EFS_mean_byenc_per_label] =  compute_stats_table(label_table)
% Step 1: Add a new column to represent the correct and incorrect trial groups
label_table = compute_correctness_group(label_table);

% Remove rows classified as "Other" if you only want to keep correct/incorrect trials
label_table = label_table(~strcmp(label_table.correctness_group, "Other"), :);

% Step 2: Use groupsummary to group by anatomical label and correctness group
EMS_EFS_mean_byenc_per_label = groupsummary(label_table, {'anatomical_label', 'correctness_group'}, 'mean', {'EMS_means', 'EFS_means'});

% Step 3: Calculate additional metrics for each group
EMS_EFS_mean_byenc_per_label.Average_EMS = mean(EMS_EFS_mean_byenc_per_label.mean_EMS_means, 2);
EMS_EFS_mean_byenc_per_label.Average_EFS = mean(EMS_EFS_mean_byenc_per_label.mean_EFS_means, 2);
EMS_EFS_mean_byenc_per_label.diff = EMS_EFS_mean_byenc_per_label.mean_EMS_means - EMS_EFS_mean_byenc_per_label.mean_EFS_means;
EMS_EFS_mean_byenc_per_label.avg_diff = mean(EMS_EFS_mean_byenc_per_label.diff, 2);

% Display the results
disp('EMS and EFS Means by Anatomy and Correctness Group:');
disp(EMS_EFS_mean_byenc_per_label);
end

function [difference_table] = compute_diff_btwn_corr_incorr(sorted_table)

% Initialize a new table to store the results
unique_labels = unique(sorted_table.anatomical_label);
n_labels = numel(unique_labels);

% Preallocate the new table
difference_table = table('Size', [n_labels, 4], ...
                         'VariableTypes', {'string', 'double', 'double', 'double'}, ...
                         'VariableNames', {'anatomical_label', 'avg_diff_Correct', 'avg_diff_Incorrect', 'avg_diff_Difference'});

% Loop through each unique anatomical label
for i = 1:n_labels
    % Get the current label
    current_label = unique_labels{i};
    
    % Find the rows corresponding to this label and each correctness group
    correct_row = strcmp(sorted_table.anatomical_label, current_label) & strcmp(sorted_table.correctness_group, "Correct_Trial");
    incorrect_row = strcmp(sorted_table.anatomical_label, current_label) & strcmp(sorted_table.correctness_group, "Incorrect_Trial");
    
    % Extract the avg_diff for correct and incorrect trials
    avg_diff_correct = sorted_table.avg_diff(correct_row);
    avg_diff_incorrect = sorted_table.avg_diff(incorrect_row);
    
    % Store the values in the new table
    difference_table.anatomical_label(i) = current_label;
    difference_table.avg_diff_Correct(i) = avg_diff_correct;
    difference_table.avg_diff_Incorrect(i) = avg_diff_incorrect;
    difference_table.avg_diff_Difference(i) = avg_diff_correct - avg_diff_incorrect;
end

% Calculate the percent change from correct
difference_table.Percent_Change_From_Correct = (difference_table.avg_diff_Difference ./ difference_table.avg_diff_Correct) * 100;

% Display the resulting table
disp(difference_table);
end
