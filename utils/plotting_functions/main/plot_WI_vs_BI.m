function [real_t_values, real_p_values, surrogate_t_values, surrogate_p_values, final_p] = plot_WI_vs_BI(params, plot_params)
    tic
    patient_id = [plot_params.patient_id];
    channel_ids_to_use = [plot_params.chan_id];
    patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, patient_id);

    for session_idx = 1:length(patient_preprocessed_data_paths)
        params.session_id = session_idx;
        patient_preprocessed_data_path = patient_preprocessed_data_paths{session_idx};
        
        plot_params.WI_BI_folder_to_save_in = sprintf("results/WI vs BI/p%s chan%s image%s enc%s sess%d", ...
                    num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
                    num2str(plot_params.image_id), num2str(plot_params.enc_id), session_idx);
    
        folder_to_save_in = plot_params.WI_BI_folder_to_save_in;
    
        warning('off', 'MATLAB:MKDIR:DirectoryExists');
        mkdir(folder_to_save_in);
        warning('on', 'MATLAB:MKDIR:DirectoryExists');
    
    
        % patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
        % image_labels = load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted
        if params.k12wm
            prefix = strsplit(patient_preprocessed_data_path,'\');
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
        
        for chan_idx = 1:length(channel_ids_to_use)
            chan_id = channel_ids_to_use(chan_idx);
            plot_params.anat = string(anat_labels.labelsanatbkedit.anatmacro1(chan_id));
            fprintf("   %s\n", plot_params.anat)
    
            % get WI, BI matrices; size: (nEtimes, nMtimes, nTrialCombinations)
            [WI, BI] = get_and_pre_process_WI_BI(plot_params.patient_id, params.comp_options, plot_params.enc_id, ...
                plot_params.image_id, plot_params.chan_id, params, plot_params.type, plot_params.only_all3_correct);
    
            % load some previous p vlaues if needed @DEPRECATING
            % selected_p_data = load(sprintf("results/selected_p_data_p%s_image%s_chan%s.mat", ...
            %     num2str(plot_params.patient_id),num2str(plot_params.image_id),num2str(plot_params.chan_id)));
    
            % plot EMS timecourse
            fig = plot_similarity_means_timecourse(WI,BI, plot_params);
            saveas(fig, sprintf("%s/WI_BI_time_courses_same_n.fig", folder_to_save_in))
    
            % plot EMS heatmap
            for same_n = [true, false]
                plot_params.same_n = same_n;
                fig = plot_similarity_means_heatmap(WI,BI, plot_params);
                if plot_params.same_n
                    saveas(fig, sprintf("%s/WI_BI_heatmaps_same_n.fig", folder_to_save_in))
                else
                    saveas(fig, sprintf("%s/WI_BI_heatmaps.fig", folder_to_save_in))
                end
            end

            %%
            % real_t_values = nan;
            % real_p_values=nan;
            % surrogate_t_values=nan;
            % surrogate_p_values = nan;
            % final_p = nan;
            % return

            %%
            % close all % close first 3 figures
    
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
    
            if ~exist(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in), 'file') || ~plot_params.skip_existing
                % perform t-test on real differences
                fprintf("    calculating real t-tests\n")
                [real_t_values, real_p_values] = calc_ttest(diff);
            
                % calculate 1000 surrogate differences and t-tests
                fprintf("    calculating surrogate differences and t-tests\n")
                [surrogate_t_values, surrogate_p_values, surrogate_diff] = calc_surrogate(WI, BI, 1000, plot_params.average_diff, plot_params.mean_out_time_dimensions);
                save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in),"real_p_values", "real_t_values", "surrogate_t_values", "surrogate_p_values", "surrogate_diff")
            else
                load(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in))
                % if ~exist("surrogate_diff","var")
                %     [surrogate_t_values, surrogate_p_values, surrogate_diff] = calc_surrogate(WI, BI, 1000, plot_params.average_diff, plot_params.mean_out_time_dimensions);
                %     save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in),"real_p_values", "real_t_values", "surrogate_t_values", "surrogate_p_values", "surrogate_diff")
                % end
            end
               
            % % plot t_values
            % fig = plot_t_values(real_value, null_values);
    
            % perform clustering and get final p or get final p from temporal
            % generalization
            if ~exist("final_p","var") || ~exist("all_real_clusters","var")
                if ~plot_params.mean_out_time_dimensions
                    % clustering p
                    [final_p, all_real_clusters, significant_clusters, sum_threshold, alpha] = cluster_analysis_and_visualization_abs_main(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values, plot_params);
                    save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in), "-append", "significant_clusters", "sum_threshold", "alpha", "final_p", "all_real_clusters")
                else
                    % temproal generalization
                    % final_p = mean(real_p_values > surrogate_p_values);
                    final_p = mean(real_t_values > surrogate_t_values);
                end
            end
    
            % plot significant clusters
            fig = plot_WI_BI_significant_clusters(BI,WI, significant_clusters, all_real_clusters, plot_params);
            saveas(fig, sprintf("%s/main_cluster_visual.fig", folder_to_save_in))
            saveas(fig, sprintf("%s/main_cluster_visual.png", folder_to_save_in))
    
            % plot differences distributions
            % if ~plot_params.mean_out_time_dimensions
            %     fig = plot_diff(diff, surrogate_diff(:), final_p, plot_params);
            %     if plot_params.average_diff
            %         saveas(fig, sprintf("%s/WI_minus_BI_distributions_trialavgs.fig", folder_to_save_in))
            %     else
            %         saveas(fig, sprintf("%s/WI_minus_BI_distributions_trial_individuals.fig", folder_to_save_in))
            %     end
            % end
    
        end
        fprintf("finished computing in %.2f seconds\n",toc);
        % close all
    end
end