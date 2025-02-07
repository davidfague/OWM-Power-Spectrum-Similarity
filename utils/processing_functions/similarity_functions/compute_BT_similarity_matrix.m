function similarity_matrix = compute_BT_similarity_matrix(mean_PS_vectors1, mean_PS_vectors2, window1_IDs, window2_IDs)
    % preallocate
    % similarity_matrix = zeros(length(window1_IDs), length(window2_IDs));
    % tic % vectorizing computation is faster and gives same results
    window1_vectors = mean_PS_vectors1(window1_IDs, :);
    window2_vectors = mean_PS_vectors2(window2_IDs, :);

    % Compute correlations and Fisher's z-transform
    corr_matrix = corr(window1_vectors', window2_vectors', 'type', 'Spearman');
    similarity_matrix = 0.5 * log((1 + corr_matrix) ./ (1 - corr_matrix));

end