%% all patients table
%%
hellbender = false;
skip_existing = false;
%%
% cd('D:\Power Spectrum Similarity')% cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
if hellbender
    addpath 'Raw Data Storage'               %#ok<UNRCH>
    addpath 'subfunctions'
else
    addpath('../Z_Raw Data Storage')      %#ok<UNRCH>    
    addpath('../subfunctions')
end

% if hellbender
    patient_IDs = [201907 201908, 201903, 201905, 201906, 201901, 201910, 201915]; %#ok<UNRCH>
% else
%     patient_IDs = [201907]; %#ok<NBRAK2>
% end


% output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc');
% output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
if hellbender
    output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc'); %#ok<UNRCH>
else 
    % output_folder = 'D:\Power Spectrum Similarity\parallel
    % output\allpatients gammamod allregions allitem enc1 correct'; 
    output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct'; %#ok<UNRCH>
end
%% main
all_patients_table = table();
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    PS_file = matfile(strcat(fullfile(output_folder, num2str(patient_ID)), '.mat'));
    ES_file = matfile(fullfile(strrep(PS_file.Properties.Source, '.mat', 'all3_ES.mat')));
    label_table = PS_file.label_table;

    if ~ismember('EMS_means', label_table.Properties.VariableNames) | ~skip_existing
        [EMS_means, EFS_means] = compute_mean_EMS_EFS(ES_file);
        label_table.EMS_means = EMS_means;
        label_table.EFS_means = EFS_means;
        clear EMS_means EFS_means
        save(PS_file.Properties.Source, 'label_table', '-append')
    end

    if ~ismember('EM_EF_diff', label_table.Properties.VariableNames) | ~skip_existing
        [EMS_means, EFS_means] = compute_mean_EMS_EFS(ES_file);
        label_table.EM_EF_diff = label_table.EMS_means - label_table.EFS_means;
        save(PS_file.Properties.Source, 'label_table', '-append')
    end

    % disp(patient_ID)
    % disp(size(label_table));
    % disp(label_table(1,:))
    % 
    % all_patients_table = [all_patients_table; label_table];
end