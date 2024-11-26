function [test_mat, control_mats] = pick_rand_test(all_avg_matrices)
    % Step 1: Generate a random index to select the test matrix
    random_index = randi(size(all_avg_matrices, 3));
    % Step 2: Extract the randomly selected matrix as test_mat
    test_mat = all_avg_matrices(:, :, random_index);
    % Step 3: Remove the selected test matrix from the original matrix to get control_mats
    control_mats = all_avg_matrices(:, :, [1:random_index-1, random_index+1:end]);
end