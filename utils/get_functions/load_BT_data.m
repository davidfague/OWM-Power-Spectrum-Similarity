function [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params, btwn_trial_type)
    % Load BT-BI data for this channel
    BT_BI_data = load(sprintf('%s\\%s\\session%d\\%s\\%s\\enc%s_image%s\\BT_%s.mat', ...
        params.output_folder, num2str(patient_id), params.session_id, comp_options{1}, btwn_trial_type, num2str(enc_id), num2str(image_id), num2str(chan_id)));
    clear BT_ES
    % BT_BI_data_std = control_EMS_matrices_std; % Rename for clarity

    % Load BT-WI data for this channel
    BT_WI_data = load(sprintf('%s\\%s\\session%d\\%s\\%s\\enc%s_image%s\\BT_%s.mat', ...
        params.output_folder, num2str(patient_id), params.session_id, comp_options{2}, btwn_trial_type, num2str(enc_id), num2str(image_id), num2str(chan_id)));
end