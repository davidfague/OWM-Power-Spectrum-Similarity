function [final_p, t_test_info, clusters_info] = WI_vs_BI(params, plotting, WI, BI)
    fprintf("processing WI vs BI for" + ...
        "   p%s session%d chan%s image%s %d-%dHz\n", ...
        num2str(params.patient_id), params.session_id, num2str(params.chan_id), ...
        num2str(params.image_id), params.freq_min, params.freq_max);
    tic

    if ~isfield(params, 'WI_BI_folder_to_save_in')
        params.patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, params.patient_id);
        params.patient_preprocessed_data_path = params.patient_preprocessed_data_paths{params.session_id};
        
        params.WI_BI_folder_to_save_in = fullfile(sprintf("results/WI vs BI/p%s chan%s image%s enc%s sess%d", ...
                    num2str(params.patient_id), num2str(params.chan_id), ...
                    num2str(params.image_id), num2str(params.enc_id), params.session_id));
        if ~isfield(params, 'anat')       
            params.anat_labels = get_anat_labels(params.patient_preprocessed_data_path, params);
            
            % params.image_labels = load(fullfile(params.patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted
            
            params.anat = string(params.anat_labels.labelsanatbkedit.anatmacro1(params.chan_id));
        end
    end

    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir(fullfile(params.WI_BI_folder_to_save_in));
    warning('on', 'MATLAB:MKDIR:DirectoryExists');

    % get WI, BI matrices; size: (nEtimes, nMtimes, nTrialCombinations)
    % [WI, BI] = get_and_pre_process_WI_BI(params);
    if isempty(WI) | isempty(BI)
        [WI, BI] = load_WI_BI_for_channel(params);
        if isempty(WI) | isempty(BI)
            warning('Failed to load WI_BI for %s', params.WI_BI_folder_to_save_in)
            final_p = nan;
            t_test_info = [];
            clusters_info = [];
            return
        end
    end

    if params.clip_inf_similarities
       [WI, BI] = clip_infs_of_z_similarities(WI,BI);
    end

    if numel(BI) == 0 || numel(WI) == 0
        fprintf("issue with this BI. Likely some some nans everywhere... skipping\n")
        final_p = nan;
        t_test_info = [];
        clusters_info = [];
        return
    end

    % plot EMS timecourse and heatmap
    if plotting
        plot_EMS(WI, BI, params, params.WI_BI_folder_to_save_in)
    end

    if plotting
        title_prefix = sprintf("p%s session%d chan%s image%s %d-%dHz\n", ...
        num2str(params.patient_id), params.session_id, num2str(params.chan_id), ...
        num2str(params.image_id), params.freq_min, params.freq_max);
        fig = plot_fft_WI_vs_BI(WI,BI, false, title_prefix, params);
        saveas(fig,sprintf("%s/FFTs.png", params.WI_BI_folder_to_save_in))
    end

    %% WI vs BI
    % temporal generalization ( mean out time dimensions
    if params.mean_out_time_dimensions
        WI = mean(WI, [1 2]);
        BI = mean(BI, [1 2]);
    end

    % calculate real WI-BI differences
    t_test_info = perform_all_t_tests(WI, BI, params);
       
    % % plot t_values
    % fig = plot_t_values(real_value, null_values);
    if plotting && ~params.without_precompute_diff
        fig = plot_differences(t_test_info, params);
        saveas(fig, sprintf("%s/WI_minus_BI.fig", params.WI_BI_folder_to_save_in))
    end

    if params.mean_out_time_dimensions
        if t_test_info.real_t_values == 0 || isnan(t_test_info.real_t_values)
            error("unexpected t value: %s", num2str(t_test_info.real_t_values))
        elseif t_test_info.real_t_values > 0 % positive
            final_p = mean(t_test_info.surrogate_t_values > t_test_info.real_t_values); % right sided
        else % negative
            final_p = mean(t_test_info.surrogate_t_values < t_test_info.real_t_values); % left sided
        end
        clusters_info = [];
        return % do not do the plotting at the bottom.
    else
        % perform cluster significance testing
        [clusters_info] = measure_cluster_significance( ...
            t_test_info.real_p_values, ...
            t_test_info.real_t_values, ...
            t_test_info.surrogate_p_values, ...
            t_test_info.surrogate_t_values);

        save(sprintf("%s/significant_cluster_data.mat", params.WI_BI_folder_to_save_in), "clusters_info")

        cluster_visualization(t_test_info, clusters_info, params)

        final_p = min([clusters_info.p_largest_negative_cluster, clusters_info.p_largest_positive_cluster]);
    end

    % plot significant clusters
    if plotting
        fig = plot_WI_BI_significant_clusters(BI,WI, params, clusters_info);
        saveas(fig, sprintf("%s/main_cluster_visual.fig", params.WI_BI_folder_to_save_in))
        saveas(fig, sprintf("%s/main_cluster_visual.png", params.WI_BI_folder_to_save_in))
    end

    % plot differences distributions
    % if ~params.mean_out_time_dimensions
    %     fig = plot_diff(diff, surrogate_diff(:), final_p, params);
    %     if params.average_diff
    %         saveas(fig, sprintf("%s/WI_minus_BI_distributions_trialavgs.fig", params.WI_BI_folder_to_save_in))
    %     else
    %         saveas(fig, sprintf("%s/WI_minus_BI_distributions_trial_individuals.fig", params.WI_BI_folder_to_save_in))
    %     end
    % end

    fprintf("finished computing in %.2f seconds\n",toc);
end