function diff = calc_diff_all_pairs(WI, BI, num_pairs)
    % Inputs:
    % WI - The within-item matrix (size: time x time x num_WI)
    % BI - The between-item matrix (size: time x time x num_BI)
    % num_pairs - Maximum number of pairs to compute (randomly selected if less than total_pairs)

    % Dimensions of the input matrices
    num_WI = size(WI, 3);
    num_BI = size(BI, 3);

    % Determine the total number of pairs
    total_pairs = num_WI * num_BI;

    % Validate num_pairs
    if nargin < 3 || isempty(num_pairs)
        num_pairs = total_pairs; % Default to all pairs if num_pairs not specified
    end
    num_pairs = min(num_pairs, total_pairs);

    % Generate random pair indices if num_pairs is less than total_pairs
    if num_pairs < total_pairs
        % rng(0); % For reproducibility (optional, remove if not needed)
        random_indices = randperm(total_pairs, num_pairs);
    else
        random_indices = 1:total_pairs; % Use all combinations if num_pairs equals total_pairs
    end

    % Preallocate the output difference matrix
    diff = zeros(size(WI, 1), size(WI, 2), num_pairs);

    % Calculate differences for the selected combinations
    for k = 1:num_pairs
        [i, j] = ind2sub([num_WI, num_BI], random_indices(k)); % Map linear index to subscripts
        diff(:, :, k) = WI(:, :, i) - BI(:, :, j);
    end
end