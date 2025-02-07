function [ES_matrix, label_table] = subset_data_by_criteria(ES_matrix, label_table, target_enc_ID, target_enc_correctness, target_anat, target_image_ids, target_enc_ids)
    % Subset data based on anatomical label, encoding ID, and encoding correctness
    
    % Get rows matching the target anatomical label
    rows_with_anat = subset_by_anat(label_table, target_anat);
    
    % Get rows matching the target encoding ID and correctness
    rows_with_enc_target = subset_by_enc(label_table, target_enc_ID, target_enc_correctness);
    
    % Get rows with the specified image and encoding IDs
    rows_with_image_at_enc = get_rows_with_image_at_enc(label_table, target_image_ids, target_enc_ids);
    
    % Combine criteria for rows to use
    rows_to_use = rows_with_anat & rows_with_enc_target & rows_with_image_at_enc;
    
    % Subset the ES_matrix and label_table
    [ES_matrix, label_table] = subset_data(ES_matrix, label_table, rows_to_use);
end

function rows = get_rows_with_image_at_enc(label_table, target_image_ids, target_enc_ids)
    % Filter rows based on target image and encoding IDs
    % Inputs:
    % - label_table: the table containing the fields 'encoding_correctness' and 'encID_to_imageID'
    % - target_image_ids: vector of image IDs to check
    % - target_enc_ids: vector of encoding IDs to check

    % Create a logical mask for encoding correctness for the selected target_enc_ids
    enc_correct_mask = ismember(label_table.encoding_correctness(:, target_enc_ids), 1);
    
    % Create a logical mask for matching image IDs in the selected target_enc_ids
    image_id_mask = ismember(label_table.encID_to_imageID(:, target_enc_ids), target_image_ids);

    % Combine both masks using logical AND across the selected target_enc_ids
    rows = any(enc_correct_mask & image_id_mask, 2);
end

function rows_with_anat = subset_by_anat(label_table, target_anat)
    % Subset rows based on anatomical label
    if isempty(target_anat) || strcmp(target_anat, 'all')
        rows_with_anat = true(size(label_table, 1), 1);
    else
        rows_with_anat = strcmp(label_table.anatomical_label, target_anat) | ismember(label_table.anatomical_label, target_anat);
    end
end

function rows_with_enc_target = subset_by_enc(label_table, target_enc_ID, target_enc_correctness)
    % Subset rows based on encoding ID and encoding correctness
    if isempty(target_enc_correctness) || length(target_enc_correctness) > 1
        rows_with_enc_target = true(size(label_table, 1), 1);
    else
        if length(target_enc_ID) > 1
            % Check if all specified encoding positions are 1
            rows_with_enc_target = all(label_table.encoding_correctness(:, target_enc_ID) == 1, 2);
        else
            rows_with_enc_target = label_table.encoding_correctness(:, target_enc_ID) == target_enc_correctness | ...
                ismember(target_enc_correctness, label_table.encoding_correctness(:, target_enc_ID));
        end
    end
end

function [ES_matrix, label_table] = subset_data(ES_matrix, label_table, rows)
    % Subset the ES_matrix and label_table based on the logical index 'rows'
    ES_matrix = ES_matrix(rows, :, :, :);
    label_table = label_table(rows, :);
end



% function [ES_matrix, label_table] = subset_data_by_criteria(ES_matrix, label_table, target_enc_ID, target_enc_correctness, target_anat)
% rows_with_anat = subset_by_anat(label_table, target_anat);
% 
% rows_with_enc_target = subset_by_enc(label_table, target_enc_ID, target_enc_correctness);
% 
% rows_with_image_at_enc = get_rows_with_image_at_enc(label_table, target_image_ids, target_enc_ids);
% 
% rows_to_use = rows_with_anat & rows_with_enc_target & rows_with_image_at_enc;
% 
% [ES_matrix, label_table] = subset_data(ES_matrix, label_table, rows_to_use);
% 
% end
% 
% function rows = get_rows_with_image_at_enc(label_table, target_image_ids, target_enc_ids)
%     % Function to filter rows based on target image and encoding IDs
%     % Inputs:
%     % - label_table: the table containing the fields 'encoding_correctness' and 'image_id'
%     % - target_image_ids: a vector of image IDs to check (e.g., 1:9, [2, 3], etc.)
%     % - target_enc_ids: a vector of encoding IDs to check (e.g., 1:3, [1, 3], etc.)
% 
%     % Create a logical mask for encoding correctness for the selected target_enc_ids
%     enc_correct_mask = ismember(label_table.encoding_correctness(:, target_enc_ids), 1);
% 
%     % Create a logical mask for matching image IDs in the selected target_enc_ids
%     image_id_mask = ismember(label_table.encID_to_imageID(:, target_enc_ids), target_image_ids);
% 
%     % Combine both masks using logical AND across the selected target_enc_ids
%     rows = any(enc_correct_mask & image_id_mask, 2);
% end
% 
% function rows_with_anat = subset_by_anat(label_table, target_anat)
%     if isempty(target_anat) | target_anat == 'all'
%         rows_with_anat = true(size(label_table, 1));
%     else
%         rows_with_anat = label_table.anatomical_label == target_anat |ismember(label_table.anatomical_label, target_anat);
%     end
% end
% 
% function rows_with_enc_target = subset_by_enc(label_table, target_enc_ID, target_enc_correctness)
%     if length(target_enc_correctness) > 1 | isempty(target_enc_correctness)
%         rows_with_enc_target = true(size(label_table),1);
%     else
%         rows_with_enc_target = label_table.encoding_correctness(:, target_enc_ID) == target_enc_correctness|ismember(target_enc_correctness, label_table.encoding_correctness(:, target_enc_ID));
%         % if target_enc_ID = 1:3 or [1,2,3] will this code ensure all 3
%         % encoding positions are 1?
%     end
% end
% 
% function [ES_matrix, label_table] = subset_data(ES_matrix, label_table, rows)
%     ES_matrix = ES_matrix(rows,:,:,:);
%     label_table = label_table(rows,:);
% end