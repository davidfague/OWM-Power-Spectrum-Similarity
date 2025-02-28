% function for getting custom params
function custom_params = get_custom_params(custom_custom_params)

    custom_params = struct();
    custom_params.k12wm = false;
    custom_params.output_folder_name = 'middle_fixation_baseline';
    custom_params.hellbender = false;

    %% override if something is passed
    if nargin >= 1 && ~isempty(custom_custom_params)
        fnames = fieldnames(custom_custom_params);
        for i = 1:numel(fnames)
            % % If the custom field exists, override the default value.
            % % (Extra fields not already in params produce a warning.)
            % if isfield(params, fnames{i})
                custom_params.(fnames{i}) = custom_custom_params.(fnames{i});
            % else
            %     warning('Unknown parameter: %s. Ignored.', fnames{i});
            % end
        end
    end

end