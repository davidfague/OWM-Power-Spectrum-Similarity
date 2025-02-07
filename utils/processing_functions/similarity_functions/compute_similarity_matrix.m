function similarity_matrix = compute_similarity_matrix(mean_PS_vectors, window1_IDs, window2_IDs)
    % preallocate
    % similarity_matrix = zeros(length(window1_IDs), length(window2_IDs));
    % tic % vectorizing computation is faster and gives same results
    window1_vectors = mean_PS_vectors(window1_IDs, :);
    window2_vectors = mean_PS_vectors(window2_IDs, :);

    % Compute correlations and Fisher's z-transform
    corr_matrix = corr(window1_vectors', window2_vectors', 'type', 'Spearman');
    similarity_matrix = 0.5 * log((1 + corr_matrix) ./ (1 - corr_matrix));
    % disp(toc)
    % similarity_matrix(1:5, 1:5)

    % % Calculate the similarity matrix
    % tic
    % for i = 1:length(window1_IDs)
    %     window1_ID = window1_IDs(i);
    %     window1_mean_power_spectrum_vector = mean_PS_vectors(window1_ID, :);
    % 
    %     for j = 1:length(window2_IDs)
    %         window2_ID = window2_IDs(j);
    %         window2_mean_power_spectrum_vector = mean_PS_vectors(window2_ID, :);
    % 
    %         % Compute correlation and Fisher z-transform
    %         r = corr(window1_mean_power_spectrum_vector', window2_mean_power_spectrum_vector', 'type', 'spearman');
    %         z = 0.5 * log((1 + r) / (1 - r));
    % 
    %         % Store similarity value
    %         similarity_matrix(i, j) = z;
    %     end
    % end
    % disp(toc)
    % similarity_matrix(1:5, 1:5)
end

% function similarity_matrix = compute_similarity_matrix(mean_PS_vectors, window1_IDs, window2_IDs)
% 
%     % Extract relevant power spectrum vectors
%     window1_mean_power_spectrum_vectors = mean_PS_vectors(window1_IDs, :, :);
%     window2_mean_power_spectrum_vectors = mean_PS_vectors(window2_IDs, :, :);
% 
%     % Preallocate similarity matrix
%     similarity_matrix = zeros(length(window1_IDs), length(window2_IDs));
% 
%     % Compute all correlations at once using vectorized operations
%     for i = 1:length(window1_IDs)
%         window1_vector = squeeze(window1_mean_power_spectrum_vectors(i, :, :));  % Extract for window1
%         for j = 1:length(window2_IDs)
%             window2_vector = squeeze(window2_mean_power_spectrum_vectors(j, :, :));  % Extract for window2
% 
%             % Compute Spearman correlation
%             r = corr(window1_vector', window2_vector', 'type', 'Spearman');
% 
%             % Apply Fisher's z-transform
%             z = 0.5 * log((1 + r) / (1 - r));
% 
%             % Store the result
%             similarity_matrix(i, j) = z;
%         end
%     end
% end


% 
% function [similarity_matrix] = compute_similarity_matrix(mean_PS_vectors, window1_IDs, window2_IDs)
%     % Get sizes
%     num_win1 = length(window1_IDs);
%     num_win2 = length(window2_IDs);
% 
%     % Preallocate temporary result array
%     temp_results = zeros(num_win1 * num_win2, 1);
% 
%     % Parallel loop over all (i, j) pairs
%     for idx = 1:num_win1 * num_win2
%         % Convert linear index to subscripts
%         [i, j] = ind2sub([num_win1, num_win2], idx);
% 
%         % Extract window IDs
%         window1_ID = window1_IDs(i);
%         window2_ID = window2_IDs(j);
% 
%         % Extract the mean power spectrum vectors
%         window1_mean_power_spectrum_vector = mean_PS_vectors(window1_ID, :);
%         window2_mean_power_spectrum_vector = mean_PS_vectors(window2_ID, :);
% 
%         % Compute correlation and Fisher z-transform
%         r = corr(window1_mean_power_spectrum_vector', window2_mean_power_spectrum_vector', 'type', 'spearman');
%         z = 0.5 * log((1 + r) / (1 - r));
% 
%         % Store the result in the temporary array
%         temp_results(idx) = z;
%     end
% 
%     % Reshape the temporary results back into the 2D similarity matrix
%     similarity_matrix = reshape(temp_results, num_win1, num_win2);
% end