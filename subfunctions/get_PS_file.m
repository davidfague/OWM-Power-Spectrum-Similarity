function PS_file = get_PS_file(output_folder, patient_ID)
    PS_file = matfile(strcat(fullfile(output_folder, num2str(patient_ID)), '.mat'));
    if ~exist(PS_file.Properties.Source, 'file')
        error("PS file does not exist: %s", PS_file.Properties.Source)
    end
end