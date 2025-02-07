function [fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs, non_recall_win_IDs, all_win_IDs] = get_window_IDs(params)        
    fixation_win_IDs = false(size(params.time));
    fixation_win_IDs(1:params.stim1_start) = true;
    fixation_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(fixation_win_IDs, params.window_start_times, params.window_end_times);
    
    enc1_win_IDs = false(size(params.time));
    enc1_win_IDs(params.stim1_start:params.stim1_end) = true;
    enc1_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(enc1_win_IDs, params.window_start_times, params.window_end_times);
    
    enc2_win_IDs = false(size(params.time));
    enc2_win_IDs(params.stim2_start:params.stim2_end) = true;
    enc2_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(enc2_win_IDs, params.window_start_times, params.window_end_times);
    
    enc3_win_IDs = false(size(params.time));
    enc3_win_IDs(params.stim3_start:params.stim3_end) = true;
    enc3_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(enc3_win_IDs, params.window_start_times, params.window_end_times);
    
    maint_win_IDs = false(size(params.time));
    maint_win_IDs(params.stim3_end : (params.stim3_end + params.maintenance_duration-1)) = true;
    maint_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(maint_win_IDs, params.window_start_times, params.window_end_times);

    non_recall_win_IDs = true(size(params.time));
    non_recall_win_IDs(params.stim3_end+params.maintenance_duration:end) = false;
    non_recall_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(non_recall_win_IDs, params.window_start_times, params.window_end_times);

    all_win_IDs = true(size(params.time));
    all_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(all_win_IDs, params.window_start_times, params.window_end_times);
end