hellbender = false;
% use_gamma_mod_chans = [true];%, false]; % [true, false] means all channels
brain_anatomies_to_process = {};
skip_existing = false;
output_folder_name = 'allpatients gammamod allregions allitem allenc baseline across trials';
%%
% cd('D:\Power Spectrum Similarity')% cd('D:\Power Spectrum Similarity')%cd('C:\Users\david\Downloads\Power Spectrum Similarity')
if hellbender
    addpath 'Raw Data Storage'               %#ok<UNRCH>
    addpath 'subfunctions'
else
    addpath('../Z_Raw Data Storage')      %#ok<UNRCH>    
    addpath('../subfunctions')
end

if hellbender
    patient_IDs = [201907 201908, 201903, 201905, 201906, 201901, 201910, 201915]; %#ok<UNRCH>
else
    patient_IDs = [201901, 201910, 201907]; %#ok<NBRAK>
end

% output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc');
% output_folder = 'D:\Power Spectrum Similarity\parallel output\allpatients gammamod allregions allitem enc1 correct';
if hellbender
    output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/', output_folder_name); %#ok<UNRCH>; %#ok<UNRCH>
else 
    % output_folder = 'D:\Power Spectrum Similarity\parallel
    % output\allpatients gammamod allregions allitem enc1 correct'; 
    % output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct'; %#ok<UNRCH>
    output_folder = fullfile('..\AA_Processed Data\', output_folder_name);
end
clear output_folder_name

% if ~exist(output_folder, 'dir')
%     mkdir(output_folder)
% end
%% main
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    % PS_file = matfile(strcat(fullfile(output_folder, num2str(patient_ID)), '.mat'));
    PS_file = get_PS_file(output_folder, num2str(patient_ID));

        % % compute correlations
    compute_all_within_trial_similarities(PS_file, skip_existing)

    % combined_save_file = combine_enc_matrices(PS_file);

end
%% functions
function [] = compute_all_within_trial_similarities(PS_file, skip_existing)


    [~, enc1_win_IDs, enc2_win_IDs, enc3_win_IDs, ~, non_selection_win_IDs, ~] = get_window_IDs(); % won't change each function call
    
    % PS_file.all_windowed_mean_PS_vectors: size = all_window_IDs, frequencies, channels x trials 
    % label_table = PS_file.label_table;

    % whole_trial_ES_matrix = compute_similarity_matrix(mean_PS_vectors, enc1_win_IDs, all_win_IDs);
    fprintf('\n')

    length_table = PS_file.label_table;
    length_table = length(length_table.patient_ID);
    mean_PS_vectors = PS_file.all_windowed_mean_PS_vectors(non_selection_win_IDs,:,:);
    %% check nans
    % nanPresence = any(isnan(mean_PS_vectors), 2);
    % figure;
    % imagesc(squeeze(nanPresence));
    % colormap([1 1 1; 0 0 0]);

    %%
    encID_order = [1,3,2]; % can get rid of, was before always combining
    all3_ES_matrix = nan(length_table, length(enc1_win_IDs), length(non_selection_win_IDs), 3);
    % save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('all3_ES.mat')));
    save_file = get_ES_file(PS_file, true, false, true);
    for encID_idx = 1:3
        encID = encID_order(encID_idx);  % This will access enc1, enc3, enc2
        % save_file = fullfile(strrep(PS_file.Properties.Source, '.mat', sprintf('_ES%d.mat', encID)));
        if exist(save_file, 'file') && skip_existing
            fprintf('%s already exists. Skipping computation.\n', save_file);
            return;
        else
            fprintf('computing %s %s\n', num2str(encID), save_file);
        end

        % Determine the appropriate window IDs based on encID
        % Assign the correct window IDs based on encID
        switch encID
            case 1
                enc_win_IDs = enc1_win_IDs;
            case 2
                enc_win_IDs = enc2_win_IDs;
            case 3
                enc_win_IDs = enc3_win_IDs;
        end

        % temp_matrix = zeros(length_table, length(enc_win_IDs), length(non_selection_win_IDs), 3);
        
        % Compute similarity for all patients for this encoding ID
        for j = 1:length_table%parfor j = 1:length_table
            all3_ES_matrix(j, :, :, encID) = compute_similarity_matrix(mean_PS_vectors(:,:,j), enc_win_IDs, non_selection_win_IDs);
        end

        % save(save_file, 'temp_matrix', '-v7.3')
    end
    save(save_file, 'all3_ES_matrix', '-v7.3')
end

write_current_script_to_destination(fullfile(output_folder, num2str(patient_ID), 'WT_ES'))

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