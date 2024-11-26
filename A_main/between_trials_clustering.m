%% compute control: correct-item-encoding matched with random non-item, but all3 correct maintenance

%% Power Spectrum Similarity
% cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
% addpath 'Raw Data Storage'               
% addpath 'subfunctions'
% cd('D:\Power Spectrum Similarity\A_main\')
addpath('../Z_Raw Data Storage')          
addpath('../subfunctions')



%% Specify large data subset to process - empty means all - regions, gamma/non-gamma channels, trials, performance, items
patient_IDs =[201907];% [201907 201908, 201903, 201905, 201906, 201901, 201910, 201915];
use_gamma_mod_chans = [true];%, false]; % [true, false] means all channels
brain_anatomies_to_process = {};
% target_encodingIDs = 1;
% target_correctness = 1;
% images_to_process = {}; %P046;

% output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc');
% output_folder = 'D:\Power Spectrum Similarity\parallel output\allpatients gammamod allregions allitem enc1 correct';
output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
% note remove pwd after moving data from pwd back to the correct folder.

%'parallel output\allpatients gammamod allregions allitem enc1 correct');
%'D:\Power Spectrum Similarity\parallel output\allpatients gammamod allregions allitem enc1 correct';
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

%% constants
% [fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs] = get_window_IDs();

%% intialize parallel pool
% if isempty(gcp('nocreate'))
%     parpool('local', 2);%25); % Open a local parallel pool if none exists
% end
combined_table = table();
%% loop through patients
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    PS_file = get_PS_file(output_folder, patient_ID);

    % % % compute correlations
    for target_enc_ids = 1%:3
        for target_image_ids = 1%:9
            test_vs_ctrl(PS_file, target_enc_ids, target_image_ids);
        end
    end
end

function test_vs_ctrl(PS_file, target_enc_ids, target_image_ids)
    ctrl_file = get_ctrl_file(PS_file, target_enc_ids, target_image_ids);
    % save(all_channels_save_file, "target_image_ids", "target_enc_ids", "test_rows", "control_rows", "label_table", "all_control_EMS_matrices_by_chan", "unique_channel_IDs", "-v7.3")
    test_file = get_ES_file(PS_file);

    label_table = PS_file.label_table;
    rows_without_nan = get_rows_without_nan(label_table);
    label_table = label_table(rows_without_nan,:);

    [ES_matrix, label_table] = subset_data_by_criteria(ES_matrix, label_table, target_enc_ID, target_enc_correctness, target_anat, target_image_ids, target_enc_ids);

    % [label_table, kept_rows] = get_cleaned_data(PS_file);% Generate an array of original row indices
    % 
    % original_indices = 1:size(PS_file.label_table, 1);
    % kept_indices = original_indices(kept_rows);% Generate an array of original row indices
    % 
    % ES_matrix = ES_file.all3_ES_matrix(kept_indices, :, :, :);  % 12180          41         640           3
    % ES_matrix = ES_matrix(kept_rows,:,:,:);
    % indices = find(mf.test_rows);
    % ES_matrix = ES_matrix(mf.test_rows);
    % % test = get_test(ES_file, label_table, mf.test_rows);

end