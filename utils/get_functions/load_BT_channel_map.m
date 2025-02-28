% for channel maps
function [channel_map] = load_BT_channel_map(params, type_to_load)
% type_to_load must be either 'WI' or 'BI'
    % Load BT-BI data for this channel
    if strcmp(type_to_load, 'BI')
        BT_BI_file= fullfile(sprintf('%s\\%s\\session%d\\%s\\%s\\%s\\enc%s_image%s\\all_channels_bt_es.mat', ...
            params.output_folder, num2str(params.patient_id), params.session_id, params.comp_options{1}, params.btwn_trial_type, ...
            sprintf("%d-%d",params.freq_min,params.freq_max), ...
            num2str(params.enc_id), num2str(params.image_id))); 
        if params.hellbender
            BT_BI_file = strrep(BT_BI_file, '\', '/'); % linux insetad of windows
        end
        channel_map = load(BT_BI_file);
    elseif strcmp(type_to_load, 'WI')
        % BT_BI_data_std = control_EMS_matrices_std; % Rename for clarity
    
        % Load BT-WI data for this channel
        BT_WI_file = fullfile(sprintf('%s\\%s\\session%d\\%s\\%s\\%s\\enc%s_image%s\\all_channels_bt_es.mat', ...
            params.output_folder, num2str(params.patient_id), params.session_id, params.comp_options{2}, params.btwn_trial_type, ...
            sprintf("%d-%d",params.freq_min,params.freq_max), ...
            num2str(params.enc_id), num2str(params.image_id)));
        if params.hellbender
            BT_WI_file = strrep(BT_WI_file, '\', '/');
        end
        channel_map = load(BT_WI_file);
    else
        error("type_to_load %s is NotImplemented for function load_BT_channel_maps")
    end
end