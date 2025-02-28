function WI_vs_BI(params, plotting)
    tic

    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir(params.WI_BI_folder_to_save_in);
    warning('on', 'MATLAB:MKDIR:DirectoryExists');

    % get WI, BI matrices; size: (nEtimes, nMtimes, nTrialCombinations)
    [WI, BI] = get_and_pre_process_WI_BI(params);

    % plot EMS timecourse and heatmap
    if plotting
        plot_EMS(WI, BI, params, params.WI_BI_folder_to_save_in)
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

    % perform cluster significance testing
    [clusters_info] = measure_cluster_significance( ...
        t_test_info.real_p_values, ...
        t_test_info.real_t_values, ...
        t_test_info.surrogate_p_values, ...
        t_test_info.surrogate_t_values);

    % perform clustering and get final p or get final p from temporal
    % generalization
    if ~exist("final_p","var") || ~exist("all_real_clusters","var")
        if ~params.mean_out_time_dimensions
            % clustering p
            [final_p, all_real_clusters, significant_clusters, sum_threshold, alpha] = cluster_analysis_and_visualization_abs_main(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values, params);
            save(sprintf("%s/significant_cluster_data.mat", params.WI_BI_params.WI_BI_folder_to_save_in), "-append", "significant_clusters", "sum_threshold", "alpha", "final_p", "all_real_clusters")
        else
            % temproal generalization
            % final_p = mean(real_p_values > surrogate_p_values);
            final_p = mean(real_t_values > surrogate_t_values);
        end
    end

    % plot significant clusters
    fig = plot_WI_BI_significant_clusters(BI,WI, significant_clusters, all_real_clusters, params);
    saveas(fig, sprintf("%s/main_cluster_visual.fig", params.WI_BI_folder_to_save_in))
    saveas(fig, sprintf("%s/main_cluster_visual.png", params.WI_BI_folder_to_save_in))

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