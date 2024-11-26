function cluster_testing(control_matrix, test_matrix)
    sanity_check(control_matrix, test_matrix)
    
    p_values = compute_p_values(control_matrix, test_matrix); % Initialize p_values with the same size as test_matrix
    
    mask = p_values < 0.05;

    cluster_analysis(mask)

end

function p_values = compute_p_values(control_matrix, test_matrix)
    p_values = zeros(size(test_matrix)); % Initialize p_values with the same size as test_matrix
    for ide = 1:size(control_matrix, 1) % encoding indices
        for idm = 1:size(control_matrix, 2) % maintenance indices
            distribution = squeeze(control_matrix(ide, idm, :)); % Extract distribution across the third dimension
            test_value = test_matrix(ide, idm); % Corresponding test matrix value
            percentile = sum(distribution > test_value) / length(distribution); % Percent of control values greater than the test value
            p_values(ide, idm) = percentile;
        end
    end
end

function cluster_analysis(mask)
    numPermutations = 1000;
    clusterSizesNull = permute_mask(numPermutations, mask);
    null_distribution = combine_nulls(clusterSizesNull);

end

function significant_clusters = cluster_size_testing(null_distribution, test_distribution)
    %% significance testing
    % null_distribution = concatenatedVector;
    % test_distribution = clusterSizes(:);
    alpha = 0.05;
    threshold = prctile(null_distribution, 100 * (1 - alpha));
    significant_clusters = test_distribution > threshold;
end

function null_distribution = combine_nulls(clusterSizesNull)
    % Initialize an empty array to hold the concatenated result
    null_distribution = [];
    % Loop through each vector in the cell array and concatenate
    for i = 1:length(clusterSizesNull)
        null_distribution = [null_distribution; clusterSizesNull{i}(:)]; % Ensure column vector
    end
end

function clusterSizesNull = permute_mask(numPermutations, mask)
    clusterSizesNull = cell(numPermutations, 1);  % Preallocate for storing cluster sizes
    % Initialize an empty array to hold the concatenated result
    for i = 1:numPermutations
        % Permute the original binary mask
        shuffledMask = mask(:);
        shuffledMask = shuffledMask(randperm(numel(shuffledMask)));
        permutedMatrix = reshape(shuffledMask, size(mask));
        
        % Identify clusters in the permuted matrix
        conn = bwconncomp(permutedMatrix);
        
        % Store the sizes of clusters in this permuted matrix
        clusterSizesNull{i} = cellfun(@numel, conn.PixelIdxList);
    end
end

function sanity_check(control_matrix, test_matrix)
    % Check that ctrl_matrix has 3 dimensions
    if ndims(control_matrix) ~= 3
        error('ctrl_matrix should have 3 dimensions, but it has %d dimensions.', ndims(control_matrix));
    end
    % Check that test_matrix has 2 dimensions
    if ndims(test_matrix) ~= 2
        error('test_matrix should have 2 dimensions, but it has %d dimensions.', ndims(test_matrix));
    end
    % Check that the first two dimensions of ctrl_matrix match the size of test_matrix
    if size(control_matrix, 1) ~= size(test_matrix, 1) || size(control_matrix, 2) ~= size(test_matrix, 2)
        error('The first two dimensions of ctrl_matrix (%d, %d) must match the size of test_matrix (%d, %d).', ...
              size(control_matrix, 1), size(control_matrix, 2), size(test_matrix, 1), size(test_matrix, 2));
    end 
end