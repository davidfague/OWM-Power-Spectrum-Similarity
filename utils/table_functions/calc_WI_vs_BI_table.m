function all_p_table = calc_WI_vs_BI_table(params, all_p_table)
    if ~params.mean_out_time_dimensions
        warning("computing all_p_table without meaning out time dimensions. This could take a while.")
    end
    
    % load patient's preprocessed data
    % patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
    patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, params.patient_id);

    for session_idx = 1:length(patient_preprocessed_data_paths)
        params.session_id = session_idx;
        patient_preprocessed_data_path = patient_preprocessed_data_paths{session_idx};
        disp(patient_preprocessed_data_path)

        anat_labels = get_anat_labels(patient_preprocessed_data_path, params);
    
        for image_id = 1:9
            params.image_id = image_id;
            % channel_ids_to_use = extractIntegersFromFilenames(params.patient_id, ...
            %     params.comp_options, params.enc_id, image_id, params); 

            [WI_channel_map] = load_BT_channel_map(params, 'WI');
            [BI_channel_map] = load_BT_channel_map(params, 'BI');
            channel_ids_to_use = get_channels_from_both_channel_maps(WI_channel_map, BI_channel_map);
    
            for chan_idx = 1:length(channel_ids_to_use)
                params.chan_id = channel_ids_to_use(chan_idx);
                params.anat = string(anat_labels.labelsanatbkedit.anatmacro1(chan_id));

                WI = WI_channel_map(params.chan_id).matrix;
                BI = BI_channel_map(params.chan_id).matrix;
                [EMS_WI, EMS_BI] = pre_process_WI_BI(params, WI, BI);
                clear WI BI
        
                [final_p, ~, ~] = WI_vs_BI(params, false, EMS_WI, EMS_BI);

                % Create the row of the table
                result_table = table(params.patient_id, ...
                                 params.chan_id, final_p, params.anat, ...
                                 params.image_id, params.session_id,...
                                 'VariableNames', {'patient_id', 'chan_id', 'p', 'anat', 'image_id', 'session_id'});
            
                all_p_table = [all_p_table; result_table]; % append the row to the table
            end
        end
    end
end