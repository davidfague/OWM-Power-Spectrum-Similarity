function sanity_check_matrices(test_mat, control_mats)
    % Check that ctrl_matrix has 3 dimensions
    if ndims(control_mats) ~= 3
        error('ctrl_matrix should have 3 dimensions, but it has %d dimensions.', ndims(control_mats));
    end
    % Check that test_matrix has 2 dimensions
    if ndims(test_mat) ~= 2
        error('test_matrix should have 2 dimensions, but it has %d dimensions.', ndims(test_mat));
    end
    % Check that the first two dimensions of ctrl_matrix match the size of test_matrix
    if size(control_mats, 1) ~= size(test_mat, 1) || size(control_mats, 2) ~= size(test_mat, 2)
        error('The first two dimensions of ctrl_matrix (%d, %d) must match the size of test_matrix (%d, %d).', ...
              size(control_mats, 1), size(control_mats, 2), size(test_mat, 1), size(test_mat, 2));
    end
end
