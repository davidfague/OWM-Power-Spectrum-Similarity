
function patient_preprocessed_data_path = get_patient_preprocessed_data_path(params, patient_ID)
    if params.k12wm
        if patient_ID < 10 % For patient_ID 1 through 9, pad with two zeros.
            patientFolder = sprintf('k12wm00%s', num2str(patient_ID));
            filePattern   = sprintf('k12wm00%s_owm_load3_s*', num2str(patient_ID));
        elseif patient_ID < 100 && patient_ID >= 10 % For patient_ID 10 through 99, pad with one zero.
            patientFolder = sprintf('k12wm0%s', num2str(patient_ID));
            filePattern   = sprintf('k12wm0%s_owm_load3_s*', num2str(patient_ID));
        elseif patient_ID >= 100 && patient_ID < 1000 % For patient_ID 100 through 999, no padding is needed.
            patientFolder = sprintf('k12wm%s', num2str(patient_ID));
            filePattern   = sprintf('k12wm%s_owm_load3_s*', num2str(patient_ID));
        else
            error("patient_ID %s must be less than 1000 for k12wm data", patient_ID)
        end
        % Build the full path to the patient directory.
        patientDir = fullfile(params.preprocessed_data_location, patientFolder);

        % List all files in the patient directory that match the file pattern.
        files = dir(fullfile(patientDir, filePattern));
        
        % Extract file names into a cell array.
        fileNames = {files.name};
        
        % Remove any file names that contain 'stim'.
        validFiles = fileNames(~contains(fileNames, 'stim') & ~contains(fileNames, 'mat'));
        
        % Create a cell array of full paths.
        patient_preprocessed_data_path = cellfun(@(f) fullfile(patientDir, f), ...
                                                 validFiles, 'UniformOutput', false);
    else
        % for utah data instead
        patient_preprocessed_data_path = { fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_ID)))};
    
    end
end
