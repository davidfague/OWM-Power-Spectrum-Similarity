% Example Script for Plot and Save RSA Analysis

% Step 1: Define Input Parameters
% Specify where the power spectrum similarity data is stored
input_folder = 'D:\Power Spectrum Similarity\parallel output\allpatients gammamod allregions allitem enc1 correct';
output_folder = fullfile(input_folder, 'Processed_RSA_Results');

% Create the output folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Step 2: Specify Patient IDs
% You can define multiple patient IDs if needed
patient_IDs = [201907]; 

% Step 3: Loop Through Each Patient and Process RSA Analysis
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    
    % Load the patient's data file
    PS_file = matfile(fullfile(input_folder, sprintf('%d.mat', patient_ID)));
    
    % Call the plot_and_save_RSA function for this patient
    plot_and_save_RSA(PS_file, output_folder);
    
    fprintf('RSA analysis and plots saved for Patient %d.\n', patient_ID);
end

function plot_and_save_RSA(PS_file, output_folder)
    % PLOT_AND_SAVE_RSA - This function filters, processes, and saves results for Power Spectrum Similarity RSA analysis
    
    % Step 1: Extract the label table and remove NaN rows
    label_table = PS_file.label_table;
    rows_with_nan = any(isnan(label_table.EMS_means), 2);
    label_table = label_table(~rows_with_nan, :);
    
    % Step 2: Subset the table based on encoding performance
    label_table = subset_table_by_enc_performance(label_table);
    
    % Step 3: Extract PSVs and plot whole-trial PSVs for each anatomy and region
    % Get all windowed mean PS vectors and filter by regions/channels/trials as necessary
    rows_without_nan = find(~rows_with_nan);
    all_PS_vectors = PS_file.all_windowed_mean_PS_vectors(:, :, :);
    all_PS_vectors = all_PS_vectors(:, :, rows_without_nan);
    
    % Step 4: Plot PSVs by anatomical regions and save plots
    anatomical_labels = unique(string(label_table.anatomical_label(:)));
    for i = 1:length(anatomical_labels)
        anat_label = anatomical_labels{i};
        
        % Filter by anatomical region
        region_indices = strcmp(label_table.anatomical_label, anat_label);
        region_PSVs = all_PS_vectors(:, :, region_indices);
        
        % Plot mean PSVs for this region
        figure;
        plot(mean(region_PSVs, 3, 'omitnan'));
        title(['Whole-Trial PSVs for ', anat_label]);
        xlabel('Time Points');
        ylabel('Power Spectrum');
        
        % Save the plot
        plot_filename = fullfile(output_folder, ['PSV_', anat_label, '.png']);
        fprintf('Finish %s', plot_filename)
        saveas(gcf, plot_filename);
        close(gcf);
    end
    
    % Step 5: Calculate EMS and EFS for the cleaned data
    [EMS, EFS] = compute_EMS_EFS(PS_file, all_PS_vectors);
    
    % Step 6: Save results to the output folder
    results_file = fullfile(output_folder, ['RSA_Results_', num2str(PS_file.patient_ID), '.mat']);
    save(results_file, 'EMS', 'EFS', 'label_table', 'anatomical_labels');
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