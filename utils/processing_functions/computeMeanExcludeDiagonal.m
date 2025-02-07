function result = computeMeanExcludeDiagonal(matrix4D)
    % Validate input
    if ndims(matrix4D) ~= 4
        error('Input must be a 4D matrix.');
    end

    % Get the size of the matrix
    [dim1, dim2, nTrials1, nTrials2] = size(matrix4D);

    % Check if the last two dimensions are equal
    if nTrials1 ~= nTrials2
        error('The last two dimensions must be of equal size.');
    end

    % Create a logical mask to exclude diagonal elements
    trialMask = ~eye(nTrials1);

    % Expand the mask to match the size of the matrix
    mask = repmat(trialMask, [1, 1, dim1, dim2]);
    mask = permute(mask, [3, 4, 1, 2]);

    % Apply the mask to exclude diagonal elements
    filteredMatrix = matrix4D(mask);

    % Compute the mean of the filtered matrix
    result = mean(filteredMatrix(:), 'omitnan');
end
