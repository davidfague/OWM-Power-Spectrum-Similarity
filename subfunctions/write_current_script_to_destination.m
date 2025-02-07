%% v3
function write_current_script_to_destination(destination_path, script_path, ask_overwrite)
    if nargin < 3
        ask_overwrite = false;
    end

    % Check if the script path is valid
    if ~isfile(script_path)
        error('The specified script does not exist: %s', script_path);
    end
    
    % Ensure destination_path is a single string or character vector
    if isstring(destination_path) && numel(destination_path) > 1
        error('destination_path must be a single string or character vector.');
    elseif iscell(destination_path) && numel(destination_path) > 1
        error('destination_path must not be a cell array with multiple elements.');
    end
    
    % Convert destination_path to a character vector if needed
    if isstring(destination_path)
        destination_path = char(destination_path);
    end

    % Extract the script name and its extension
    [~, script_name, script_ext] = fileparts(script_path);
    
    % Get the current time
    currentTime = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm__');
    currentTimeStr = char(currentTime); % Convert to a character vector for use in filenames

    % Append the script name and extension to the destination path
    destination_file = fullfile(destination_path, strcat(currentTimeStr, script_name, script_ext));
    
    % Check if the destination file already exists
    if isfile(destination_file) && ask_overwrite
        prompt = sprintf('File %s already exists. Overwrite? (y/n): ', destination_file);
        user_response = input(prompt, 's');
        if ~strcmpi(user_response, 'y')
            fprintf('Operation canceled. File not overwritten.\n');
            return; % Exit the function without writing
        end
    end

    % Read the content of the script
    file_content = fileread(script_path);
    
    % Write the content to the destination file
    fid = fopen(destination_file, 'w');
    if fid == -1
        error('Failed to open the destination file for writing.');
    end
    
    % Write the content and close the file
    fwrite(fid, file_content);
    fclose(fid);
    
    fprintf('Script successfully written to: %s\n', destination_file);
end

%% v2
% function write_current_script_to_destination(destination_path)
%     % Get the name of the currently running script
%     current_file = matlab.desktop.editor.getActiveFilename();
% 
%     % Check if a file is open
%     if isempty(current_file)
%         error('No active script is currently open in the MATLAB editor.');
%     end
% 
%     % Extract the script name and its extension
%     [~, script_name, script_ext] = fileparts(current_file);
% 
%     % Append the script name and extension to the destination path
%     destination_file = fullfile(destination_path, [script_name, script_ext]);
% 
%     % Read the content of the current script
%     file_content = fileread(current_file);
% 
%     % Write the content to the destination file
%     fid = fopen(destination_file, 'w');
%     if fid == -1
%         error('Failed to open the destination file for writing.');
%     end
% 
%     % Write the content and close the file
%     fwrite(fid, file_content);
%     fclose(fid);
% 
%     fprintf('Current script successfully written to: %s\n', destination_file);
% end

%% v1
% function write_current_script_to_destination(destination_path)
%     % Get the name of the currently running script
%     current_file = matlab.desktop.editor.getActiveFilename();
% 
%     % Check if a file is open
%     if isempty(current_file)
%         error('No active script is currently open in the MATLAB editor.');
%     end
% 
%     % Read the content of the current script
%     file_content = fileread(current_file);
% 
%     % Write the content to the destination file
%     fid = fopen(destination_path, 'w');
%     if fid == -1
%         error('Failed to open the destination file for writing.');
%     end
% 
%     % Write the content and close the file
%     fwrite(fid, file_content);
%     fclose(fid);
% 
%     fprintf('Current script successfully written to: %s\n', destination_path);
% end