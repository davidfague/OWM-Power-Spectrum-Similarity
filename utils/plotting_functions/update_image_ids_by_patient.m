function updatedTable = update_image_ids_by_patient(tableData, patientColumnName, imageColumnName, numImagesPerPatient)
    % Set default values for column names and num_images_per_patient
    if nargin < 2
        patientColumnName = 'patient_id';
    end
    if nargin < 3
        imageColumnName = 'image_id';
    end
    if nargin < 4
        numImagesPerPatient = 9;
    end

    % Get unique patient IDs
    uniquePatients = unique(tableData.(patientColumnName));
    
    % Iterate over each patient and update the image_id
    for i = 1:length(uniquePatients)
        % Get the rows for the current patient
        patientRows = tableData.(patientColumnName) == uniquePatients(i);
        
        % Update the image_id for the current patient
        tableData.(imageColumnName)(patientRows) = tableData.(imageColumnName)(patientRows) + (i - 1) * numImagesPerPatient;
    end

    % Return the updated table
    updatedTable = tableData;
end
