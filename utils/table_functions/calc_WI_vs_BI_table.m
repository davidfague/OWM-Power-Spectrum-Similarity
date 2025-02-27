function all_p_table = calc_WI_vs_BI_table(params, plot_params, all_p_table)
    patient_id = [plot_params.patient_id];

        % load patient's preprocessed data
    % patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
    patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, patient_id);

    for session_idx = 1:length(patient_preprocessed_data_paths)
        params.session_id = session_idx;
        patient_preprocessed_data_path = patient_preprocessed_data_paths{session_idx};
        disp(patient_preprocessed_data_path)

        if params.k12wm
            if params.hellbender
                prefix = strsplit(patient_preprocessed_data_path,'/'); % linux
            else
                prefix = strsplit(patient_preprocessed_data_path,'\');
            end
            prefix = prefix{end};
            labels = load(fullfile(patient_preprocessed_data_path, sprintf("%s_labelsAnat.mat", prefix)));
            anat_labels = struct();
            anat_labels.labelsanatbkedit = labels.bipolarAnat;
            clear labels prefix
        else
            % image_labels = load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted
            anat_labels = load(fullfile(patient_preprocessed_data_path, ...
            "D_OWM_t_bipolar.mat"), 'labelsanatbkedit'); % i've already made sure this is the same as the table
        end
    
        for image_id = 1:9
            plot_params.image_id = image_id;
            channel_ids_to_use = extractIntegersFromFilenames(patient_id, params.comp_options, plot_params.enc_id, image_id, params); % image_id=1 is used as a default here. Channels should change with respect to iamge_id.
    
            for chan_idx = 1:length(channel_ids_to_use)
                chan_id = channel_ids_to_use(chan_idx);
                plot_param.chan_id = chan_id;
                fprintf("p%s session%d chan%s image%s %d-%dHz\n", num2str(patient_id), session_idx, num2str(chan_id), num2str(image_id), params.freq_min, params.freq_max);
                anat = string(anat_labels.labelsanatbkedit.anatmacro1(chan_id));

                % get WI, BI matrices; size: (nEtimes, nMtimes, nTrialCombinations)
                [WI, BI] = get_and_pre_process_WI_BI(patient_id, params.comp_options, plot_params.enc_id, ...
                    image_id, chan_id, params, plot_params.type, plot_params.only_all3_correct);
    
                if numel(BI) == 0 || numel(WI) == 0
                    fprintf("issue with this BI. Likely some some nans everywhere... skipping\n")
                    final_p = nan;
                    result_table = table(patient_id, ...
                                     chan_id, final_p, anat, ...
                                     image_id, params.session_id,...
                                     'VariableNames', {'patient_id', 'chan_id', 'p', 'anat', 'image_id', 'session_id'});
                
                    all_p_table = [all_p_table; result_table];
                    continue
                end
        
                if params.clip_inf_similarities
                   [WI, BI] = clip_infs_of_z_similarities(WI,BI);
                end
        
                % temporal generalization
                if plot_params.mean_out_time_dimensions
                    WI=mean(WI, [1 2]);
                    BI=mean(BI, [1 2]);
                end
        
                % calculate real differences
                fprintf("    calculating real differences\n")
                if plot_params.average_diff
                    diff = calc_diff_avg(WI, BI, 1000);
                else
                    diff = calc_diff_all_pairs(WI, BI, 10000);
                end
        
                fprintf("    calculating real t-tests\n")
                [real_t_values, real_p_values] = calc_ttest(diff);
            
                % calculate 1000 surrogate differences and t-tests
                fprintf("    calculating surrogate differences and t-tests\n")
                [surrogate_t_values, surrogate_p_values, surrogate_diff] = calc_surrogate(WI, BI, 1000, plot_params.average_diff, plot_params.mean_out_time_dimensions);
                % save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in),"real_p_values", "real_t_values", "surrogate_t_values", "surrogate_p_values", "surrogate_diff")
        
                if any(isnan(real_t_values))
                    error("nans detected in real_t_values.")
                elseif any(isnan(real_p_values))
                    error("nans detecte in real_p_values")
                elseif any(isnan(surrogate_p_values))
                    error("nans detected in surrogate_p_values")
                elseif any(isnan(surrogate_t_values))
                    error("nans detected in surrogate_t_values")
                end

                final_p = mean(real_t_values > surrogate_t_values);
        
                    % Create the table
                result_table = table(patient_id, ...
                                     chan_id, final_p, anat, ...
                                     image_id, ...
                                     session_idx, ...
                                     'VariableNames', {'patient_id', 'chan_id', 'p', 'anat', 'image_id', 'session_id'});
                
                all_p_table = [all_p_table; result_table];
            end
        end
    end
end