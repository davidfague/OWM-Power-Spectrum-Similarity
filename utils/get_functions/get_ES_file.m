% function ES_file = get_ES_file(PS_file)
%     ES_file = matfile(fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('all3_ES.mat'))));
%     if ~exist(ES_file.Properties.Source, 'file')
%         ES_file = matfile(fullfile(strrep(PS_file.Properties.Source, '.mat', 'all3Enc_ES.mat')));
%     end
% 
%     if ~exist(ES_file.Properties.Source, 'file')
%         error("ES file does not exist: %s", ES_file.Properties.Source)
%     end
% end
%%

function ES_file = get_ES_file(PS_file, input_matfile, return_matfile, ignore_exists)
    % Check if input_matfile or return_matfile is provided, default to true
    if nargin < 2
        input_matfile = true;
    end
    if nargin < 3
        return_matfile = true;
    end
    if nargin < 4
        ignore_exists = false;
    end

    % If input_matfile is true, PS_file is expected to be a matfile object
    if input_matfile
        if ~(class(PS_file) == "matlab.io.MatFile")
            error("PS_file must be a matfile when input_matfile is true.");
        end
        % Get the source file path from the matfile object
        file_path = strrep(PS_file.Properties.Source, 'PSVs.mat', 'encoding_similarity.mat');
    else
        % If input_matfile is false, PS_file is assumed to be a file path string
        if ~ischar(PS_file) && ~isstring(PS_file)
            error("PS_file must be a file path when input_matfile is false.");
        end
        % file_path = strrep(PS_file, 'PSVs.mat', sprintf('%sall3_ES.mat',num2str(patient_ID)));
        file_path = strrep(PS_file, 'PSVs.mat', 'encoding_similarity.mat');
    end

    % Check if the ES file exists
    if ~exist(file_path, 'file') && ~ignore_exists
        error("ES file does not exist: %s", file_path);
    end

    % Return either the matfile or the file path
    if return_matfile
        ES_file = matfile(file_path);
    else
        ES_file = file_path;
    end
end
