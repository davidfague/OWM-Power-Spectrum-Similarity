
%% v2 implement get_parameters.m
close all
clear all

custom_params = struct();
custom_params.k12wm = false;
custom_params.patient_IDs = [201908, 201910];

params = get_parameters(custom_params);

params.skip_existing = true; % temp param

%% main
for idx = 1:length(params.patient_IDs)
    patient_ID = params.patient_IDs(idx);
    % PS_file = matfile(strcat(fullfile(output_folder, num2str(patient_ID)), '.mat'));
    PS_file = get_PS_file(params, num2str(patient_ID));

        % % compute correlations
    compute_all_within_trial_similarities(PS_file, params)

    write_current_script_to_destination(fullfile(params.output_folder, num2str(patient_ID)), strcat(mfilename('fullpath'), '.m'));

end

%% functions
function [] = compute_all_within_trial_similarities(PS_file, params)
    
    % PS_file.all_windowed_mean_PS_vectors: size = all_window_IDs, frequencies, channels x trials 
    % label_table = PS_file.label_table;

    % whole_trial_ES_matrix = compute_similarity_matrix(mean_PS_vectors, enc1_win_IDs, all_win_IDs);
    fprintf('\n')

    length_table = PS_file.label_table;
    length_table = length(length_table.patient_ID);
    mean_PS_vectors = PS_file.all_windowed_mean_PS_vectors(params.non_recall_win_IDs,:,:);
    %% check nans
    % nanPresence = any(isnan(mean_PS_vectors), 2);
    % figure;
    % imagesc(squeeze(nanPresence));
    % colormap([1 1 1; 0 0 0]);

    %%
    all3_ES_matrix = nan(length_table, length(params.enc1_win_IDs), length(params.non_recall_win_IDs), 3);
    % save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('all3_ES.mat')));
    save_file = get_ES_file(PS_file, true, false, true);
    for encID = 1:3
        % save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('_ES%d.mat', encID)));
        if exist(save_file, 'file') && params.skip_existing
            fprintf('%s already exists. Skipping computation.\n', save_file);
            return;
        else
            fprintf('computing enc%s %s\n', num2str(encID), save_file);
        end

        % Determine the appropriate window IDs based on encID
        % Assign the correct window IDs based on encID
        switch encID
            case 1
                enc_win_IDs = params.enc1_win_IDs;
            case 2
                enc_win_IDs = params.enc2_win_IDs;
            case 3
                enc_win_IDs = params.enc3_win_IDs;
        end

        % temp_matrix = zeros(length_table, length(enc_win_IDs), length(params.non_recall_win_IDs), 3);
        
        % Compute similarity for all patients for this encoding ID
        for j = 1:length_table%parfor j = 1:length_table
            all3_ES_matrix(j, :, :, encID) = compute_similarity_matrix(mean_PS_vectors(:,:,j), enc_win_IDs, params.non_recall_win_IDs);
        end

        % save(save_file, 'temp_matrix', '-v7.3')
    end
    save(save_file, 'all3_ES_matrix', '-v7.3')
end


%% old code for computing individual and then combining; now just combine every time.
% function [combined_save_file] = combine_enc_matrices(PS_file)
% % can probably update to delete the 3 files that got combined
%     % Precompute window IDs (won't change for each function call)
%     % [fixation_win_IDs, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, maint_win_IDs, non_selection_win_IDs, all_win_IDs] = get_window_IDs();
% 
%     % Output save file for the combined matrix
%     combined_save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', '/all3 ES/all3_ES.mat'));
% 
%     if exist(combined_save_file, 'file')
%         fprintf('%s already exists. Skipping computation.\n', combined_save_file);
%         return;
%     end
% 
%     % Prepare a matrix to hold the combined similarities for each encoding period
%     all3_ES_matrix = [];
% 
%     for encID = 1:3
%         save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('_ES%d.mat', encID)));
% 
%         if exist(save_file, 'file')
%             fprintf('Loading %s\n', save_file);
%             load(save_file, 'temp_matrix');
%         else
%             % fprintf('ES%d file not found. Skipping this encoding.\n', encID);
%             % continue;
%             error('ES%d file not found. Skipping this encoding.\n', encID);
%         end
% 
%         % Add a new dimension for each encoding matrix
%         if isempty(all3_ES_matrix)
%             % Initialize the combined matrix with the first encoding
%             all3_ES_matrix = zeros(size(temp_matrix, 1), size(temp_matrix, 2), size(temp_matrix, 3), 3);
%         end
% 
%         % Store the temp_matrix in the appropriate encoding ID slice
%         all3_ES_matrix(:,:,:,encID) = temp_matrix;
%     end
% 
%     % Save the combined matrix into a new file
%     fprintf('Saving combined encoding matrices to %s\n', combined_save_file);
%     save(combined_save_file, 'all3_ES_matrix', '-v7.3');
% 
% end

%%
    % If the matrix already exists, skip execution
    % some old code for computing combined all encoding at once (prefer not to keep
    % so much in memory while computing so save 1 at a time.)
    % if isprop(save_file, 'all_all3Enc_wholeTrial_ES_matrix')
    %     disp('all_all3Enc_wholeTrial_ES_matrix already exists. Skipping computation.');
    %     return;
    % else
    %     fprintf('computing %s all_all3Enc_wholeTrial_ES_matrix', string(PS_file));
    % end
            %all_all3Enc_wholeTrial_ES_matrix = zeros(length(label_table.patient_ID) ,size(whole_trial_ES_matrix, 1),size(whole_trial_ES_matrix, 2), 3);
        %        parfor i=1:length(label_table.patient_ID) % every channel-trial combination    
        %     % Temporary variable to hold the similarity matrices for this iteration
        %     all3Enc_wholeTrial_ES_matrix = zeros(size(whole_trial_ES_matrix, 1), size(whole_trial_ES_matrix, 2), 3);
        %     % Compute similarity matrices for each encoding period
        %     all3Enc_wholeTrial_ES_matrix(:, :, 1) = compute_similarity_matrix(mean_PS_vectors, enc1_win_IDs, all_win_IDs);
        %     all3Enc_wholeTrial_ES_matrix(:, :, 2) = compute_similarity_matrix(mean_PS_vectors, enc2_win_IDs, all_win_IDs);
        %     all3Enc_wholeTrial_ES_matrix(:, :, 3) = compute_similarity_matrix(mean_PS_vectors, enc3_win_IDs, all_win_IDs);
        % 
        %     % all3Enc_EMS_matrix = all3Enc_wholeTrial_ES_matrix(:,maint_win_IDs,:);
        %     % all3Enc_EFS_matrix = all3Enc_wholeTrial_ES_matrix(:,fixation_win_IDs,:);
        % 
        %     % all3Enc_EES_matrix = all3Enc_wholeTrial_ES_matrix(:,fixation_win_IDs,:);
        % 
        %     % Assign the temporary result to the main matrix
        %     all_all3Enc_wholeTrial_ES_matrix(i, :, :, :) = all3Enc_wholeTrial_ES_matrix;
        % end