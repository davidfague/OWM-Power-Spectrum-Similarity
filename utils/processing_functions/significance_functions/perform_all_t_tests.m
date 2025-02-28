function t_test_info = perform_all_t_tests(WI, BI, params)
    file_to_save = fullfile(sprintf('%s/significant_cluster_data.mat', params.WI_BI_folder_to_save_in));

    file_exists = exist(file_to_save, 'file');
    if ~file_exists
        run_t_test = true;
    else
        run_t_test = ~any(strcmp('t_test_info', who('-file', file_to_save))) || ~params.skip_existing;
    end

    if run_t_test
        % calculate real WI-BI differences
        t_test_info = struct();
        fprintf("    calculating real WI-BI differences\n")
        if params.average_diff
            t_test_info.real_diff = calc_diff_avg(WI, BI, 1000);
        else
            t_test_info.real_diff = calc_diff_all_pairs(WI, BI, 10000);
        end

        if any(isnan(t_test_info.real_diff), "all")
            error("nan in real_diff for t-test")
        end
    
        % perform t-test on real differences
        fprintf("    calculating real t-tests\n")
        [t_test_info.real_t_values, t_test_info.real_p_values] = calc_ttest(t_test_info.real_diff);
    
        % calculate 1000 surrogate differences and t-tests
        fprintf("    calculating surrogate differences and t-tests\n")
        [t_test_info.surrogate_t_values, t_test_info.surrogate_p_values, t_test_info.surrogate_diffs] = calc_surrogate(WI, BI, 1000, params.average_diff, params.mean_out_time_dimensions);

        save(fullfile(params.WI_BI_folder_to_save_in, 'significant_cluster_data.mat'), 't_test_info')
    else
        fprintf("    skipping t_test calculation\n    loading %s", file_to_save)
        load(file_to_save, 't_test_info')
    end
end