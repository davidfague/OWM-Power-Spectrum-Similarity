function boolean_array = bool_mask_array(array_to_bool, max_index)
    % note array_to_bool should be like [3, 5, 8]; max_index = 10
    % result would be [0, 0, 1, 0, 1, 0, 0, 1, 0, 0]
    
    boolean_array = false(1, max_index);  % Initialize with false
    boolean_array(array_to_bool) = true;       % Set true at indices in sigchans
end