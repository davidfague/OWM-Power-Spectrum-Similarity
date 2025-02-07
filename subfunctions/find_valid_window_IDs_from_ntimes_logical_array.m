function valid_window_IDs = find_valid_window_IDs_from_ntimes_logical_array(ntimes_logical_array, window_start_times, window_end_times)
    num_windows = length(window_start_times);
    valid_windows = false(1, num_windows);

    for i = 1:num_windows
        if all(ntimes_logical_array(window_start_times(i):window_end_times(i)))
            valid_windows(i) = true;
        end
    end

    % Filter window1s to only use valid windows
    valid_window_IDs = find(valid_windows); 
    % num_valid_windows = length(window1_IDs);
end