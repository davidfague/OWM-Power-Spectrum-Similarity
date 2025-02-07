function [EMS_means, EFS_means] = compute_mean_EMS_EFS(ES_file)
 [fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs, non_selection_win_IDs, all_win_IDs] = get_window_IDs(); % won't change each function call
    EMS = get_enc_similarity(ES_file, maint_win_IDs);
    EFS = get_enc_similarity(ES_file, fixation_win_IDs);
    EMS_means = squeeze(compute_mean_ES(EMS));
    EFS_means = squeeze(compute_mean_ES(EFS));
end


function similarity_subset = get_enc_similarity(PS_file, select_window_ids)
    similarity_subset = zeros(size(PS_file.all3_ES_matrix, 1), size(PS_file.all3_ES_matrix, 2), length(select_window_ids), 3);
    for i = 1:size(PS_file.all3_ES_matrix, 1)
        similarity_subset(i, :,:,:) = PS_file.all3_ES_matrix(i, :, select_window_ids, :); % each encoding vs part of trial
    end
end

function means = compute_mean_ES(ES_matrix)
    % ES_means = zeros(size(ES_matrix, 1), 3);
    % for i=1:size(ES_matrix, 1) % channel-trial combinations
    mean_over_time = mean(ES_matrix, 3);

    means = mean(mean_over_time, 2); % mean across select_window_IDs then mean across encoding

    % means should be channel-trial x 3encodings
end