function [fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs, non_selection_win_IDs, all_win_IDs] = get_window_IDs()
    time = 1:9001;
    % 100ms time window with step of 10ms
    wt = length(time);
    num_windows = floor((wt - 100) / 10) + 1; % Correct calculation for number of windows
    % Precompute window start and end times
    window_start_times = time(1:10:(num_windows - 1) * 10 + 1); % Start times
    clear num_windows wt
    window_end_times = time(window_start_times + 99);            % End times
    
    % Define the time windows for each encoding period
    stim1_start = 1000; stim1_end = 1500; % stim 1 time window
    stim2_start = 1500; stim2_end = 2000; % stim 2 time window
    stim3_start = 2000; stim3_end = 2500; % stim 3 time window
    
    fixation_win_IDs = false(size(time));
    fixation_win_IDs(1:stim1_start) = true;
    fixation_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(fixation_win_IDs, window_start_times, window_end_times);
    
    enc1_win_IDs = false(size(time));
    enc1_win_IDs(stim1_start:stim1_end) = true;
    enc1_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(enc1_win_IDs, window_start_times, window_end_times);
    
    enc2_win_IDs = false(size(time));
    enc2_win_IDs(stim2_start:stim2_end) = true;
    enc2_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(enc2_win_IDs, window_start_times, window_end_times);
    
    enc3_win_IDs = false(size(time));
    enc3_win_IDs(stim3_start:stim3_end) = true;
    enc3_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(enc3_win_IDs, window_start_times, window_end_times);
    
    maint_win_IDs = false(size(time));
    maint_win_IDs(stim3_end:stim3_end+4000-1) = true;
    maint_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(maint_win_IDs, window_start_times, window_end_times);

    all_win_IDs = true(size(time));
    all_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(all_win_IDs, window_start_times, window_end_times);

    non_selection_win_IDs = true(size(time));
    non_selection_win_IDs(stim3_end+4000:end) = false;
    non_selection_win_IDs = find_valid_window_IDs_from_ntimes_logical_array(non_selection_win_IDs, window_start_times, window_end_times);

end