function ctrl_file = get_ctrl_file(PS_file, patient_id, target_enc_ids, target_image_ids)
    save_folder = strrep(PS_file.Properties.Source, '.mat', sprintf('/corrWOitemMaint vs corrcWitemEnc/between_trials_enc%d_image%d', target_enc_ids, target_image_ids));
    ctrl_file = matfile(fullfile(save_folder, 'all_channels_data.mat'));

    if ~exist(ctrl_file.Properties.Source, 'file')
        error("ctrl file does not exist %s", ctrl_file.Properties.Source)
    end
end