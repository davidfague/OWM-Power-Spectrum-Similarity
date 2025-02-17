function integersList = extractIntegersFromFilenames(patient_id, comp_options, enc_id, image_id, params)
    file_to_load = fullfile(sprintf('%s\\%s\\session%d\\%s\\%s\\%d-%d\\enc%s_image%s\\', ...
        params.output_folder, num2str(patient_id), params.session_id, comp_options{1}, params.btwn_trial_type,  params.freq_min, params.freq_max, num2str(enc_id), num2str(image_id)));
    if params.hellbender
        file_to_load = strrep(file_to_load, '\', '/'); % linux instead of windows. I thought fullfile() was supposed to automatically handle that but apparently not.
    end
    
    % Function to extract the integer at the end of each filename in a directory
    % Inputs:
    %   directoryPath - the path to the directory containing the files
    % Outputs:
    %   integersList - a list of integers extracted from the filenames

    % Get a list of all files in the directory
    files = dir(file_to_load);
    
    % Initialize the list to store integers
    integersList = [];

    % Loop through each file
    for i = 1:length(files)
        % Skip directories
        if files(i).isdir
            continue;
        end

        % Get the filename
        filename = files(i).name;

        % Split the filename by '_'
        parts = split(filename, '_');

        % Extract the last part
        lastPart = parts{end};

        % Remove file extension if present (e.g., .txt, .mat)
        [~, lastPart, ~] = fileparts(lastPart);

        % Convert to integer
        integerValue = str2double(lastPart);

        % Check if conversion is successful
        if ~isnan(integerValue)
            integersList = [integersList; integerValue];
        end
    end
end