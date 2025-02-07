function [all_t_values, all_p_values, all_diffs] = calc_surrogate(WI_data, BI_data, n_perm, average_diff, return_all_diffs)
    if nargin < 5
        return_all_diffs=false;
    end
    sanity_check_sizes(WI_data, BI_data)

    if return_all_diffs
        if average_diff
            all_diffs = nan([size(WI_data, [1 2]), 1000, n_perm]);
        else
            all_diffs = nan([size(WI_data, [1 2]), 10000, n_perm]);
        end
    else
        all_diffs = [];
    end
    
    all_t_values = nan(size(WI_data, 1), size(WI_data, 2), n_perm);
    all_p_values = nan(size(WI_data, 1), size(WI_data, 2), n_perm);
    for perm = 1:n_perm
        if perm == n_perm
            fprintf("   perm %s \n", num2str(perm))
        elseif mod(perm, 100) == 0
            fprintf("   perm %s", num2str(perm))
        end
        % shuffle WI, BI labels and calculate ttest nperm times
        [new_WI_data, new_BI_data] = shuffle_labels(WI_data, BI_data);

        % calculate differences
        if average_diff
            diff = calc_diff_avg(new_WI_data, new_BI_data, 1000);
        else
            diff = calc_diff_all_pairs(new_WI_data, new_BI_data, 10000);
        end
        if return_all_diffs
           all_diffs(:,:,:, perm) = diff;
        end

        % run t test between diff and 0
        [t_values, p_values] = calc_ttest(diff);
        all_t_values(:,:,perm) = t_values;
        all_p_values(:,:,perm) = p_values;
    end
end

function sanity_check_sizes(WI_data, BI_data)
    % make sure that first and second dimensions have same size
    assert(isequal(size(WI_data, 1), size(BI_data, 1)), 'First dimension of inputs must match.');
    assert(isequal(size(WI_data, 2), size(BI_data, 2)), 'Second dimension of inputs must match.');
end

function [new_WI_data, new_BI_data] = shuffle_labels(WI_data, BI_data)
    % first and second dimensions of inputs should be t he same.
    % shuffle slices along the 3rd dimension and place them in new
    % matrices
    combined_data = cat(3, WI_data, BI_data); % Concatenate WI_data and BI_data along the 3rd dimension (TRIAL PAIRS)

    % Shuffle slices along the 3rd dimension
    shuffled_indices = randperm(size(combined_data, 3));
    shuffled_data = combined_data(:, :, shuffled_indices);

    % Split shuffled_data back into new_WI_data and new_BI_data
    nWI = size(WI_data, 3); % Number of slices in WI_data
    new_WI_data = shuffled_data(:, :, 1:nWI);
    new_BI_data = shuffled_data(:, :, nWI+1:end);
end