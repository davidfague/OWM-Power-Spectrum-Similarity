% function PS_file = get_PS_file(output_folder, patient_ID)
%     % PS_file = matfile(strcat(fullfile(output_folder, num2str(patient_ID)), '.mat'));
%     PS_file = matfile(fullfile(strcat(fullfile(output_folder, fprintf('%s/%s.mat', num2str(patient_ID)), '.mat'))));
%     if ~exist(PS_file.Properties.Source, 'file')
%         error("PS file does not exist: %s", PS_file.Properties.Source)
%     end
% end

function PS_file = get_PS_file(params, patient_ID, use_matfile) % use_matfile returns matfile
    % Default use_matfile to true if not provided
    if nargin < 3
        use_matfile = true;
    end
    
    % Construct the file path
    % file_path = fullfile(output_folder, sprintf('%s/%s.mat', num2str(patient_ID), num2str(patient_ID)));
    if params.k12wm
        file_path = fullfile(params.output_folder, sprintf('%s/PSVs_%s_%s.mat', num2str(patient_ID), params.band_to_process.name, num2str(params.session_id)));
    else
        file_path = fullfile(params.output_folder, sprintf('%s/PSVs_%s.mat', num2str(patient_ID), params.band_to_process.name));
    end
    
    % return matfile or path
    if use_matfile
        % Check if the file exists
        if ~exist(file_path, 'file')
            error("PS file does not exist: %s", file_path);
        end
        PS_file = matfile(file_path);
    else
        % PS_file = load(file_path);
        PS_file = file_path;
    end
end
