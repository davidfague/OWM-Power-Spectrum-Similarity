
function session_ids = get_available_session_ids(params, patient_ID)
% GET_AVAILABLE_SESSION_IDS returns the session IDs available for a given patient.
%
%   session_ids = GET_AVAILABLE_SESSION_IDS(params, patient_ID) searches the folder
%   specified by params.output_folder for a subfolder corresponding to patient_ID and
%   finds all files that match the pattern 'PSVs*.mat'. It then extracts the session IDs
%   (the part after 'PSVs' and before '.mat') and returns them as a numeric array.
%
%   Inputs:
%       params     - A structure that must include at least the field:
%                       output_folder : path to the folder containing patient folders.
%       patient_ID - A number or string identifying the patient.
%
%   Output:
%       session_ids - A numeric array of session IDs available for the patient.
%
%   Example:
%       params.output_folder = '/path/to/output';
%       patient_ID = 101;
%       session_ids = get_available_session_ids(params, patient_ID);
%
%   See also: DIR, FULLFILE, STR2DOUBLE.

    % Convert patient_ID to string in case it is numeric.
    patient_folder = fullfile(params.output_folder, num2str(patient_ID));
    
    % Check if the patient folder exists.
    if ~exist(patient_folder, 'dir')
        error('Patient folder "%s" does not exist.', patient_folder);
    end

    % Build the file pattern. We expect files named like: 'PSVs<session_id>.mat'
    file_pattern = fullfile(patient_folder, 'PSVs*.mat');
    file_list = dir(file_pattern);
    
    if isempty(file_list)
        warning('No session files found for patient %s in folder %s.', num2str(patient_ID), patient_folder);
        session_ids = [];
        return;
    end
    
    % Loop through the files and extract session IDs.
    session_ids = [];  % initialize as empty
    for k = 1:length(file_list)
        filename = file_list(k).name;
        % We expect the filename to start with 'PSVs' and end with '.mat'
        prefix = 'PSVs';
        suffix = '.mat';
        if startsWith(filename, prefix) && endsWith(filename, suffix)
            % Extract the portion between the prefix and suffix.
            session_str = filename(length(prefix)+1 : end-length(suffix));
            session_num = str2double(session_str);
            if ~isnan(session_num)
                session_ids(end+1) = session_num;  %#ok<AGROW>
            else
                warning('Could not convert session ID "%s" from file "%s" to a number.', session_str, filename);
            end
        else
            warning('File "%s" does not match the expected pattern "PSVs<session_id>.mat".', filename);
        end
    end
    
    % Sort session IDs in ascending order.
    session_ids = sort(session_ids);
end