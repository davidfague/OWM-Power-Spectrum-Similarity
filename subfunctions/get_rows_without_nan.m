function rows_without_nan = get_rows_without_nan(label_table)
    % Identify rows with NaNs in the label table
    rows_without_nan = ~any(isnan(label_table.EMS_means), 2);
end

% function [label_table, all_window kept_rows] = get_cleaned_data(data_file)
%     label_table = data_file.label_table;
%     % Identify rows with NaNs in the label table
%     rows_with_nan = any(isnan(label_table.EMS_means), 2);
%     rows_without_nan = ~rows_with_nan;
%     % Filter the label table and the windowed mean PS vectors
%     label_table = label_table(rows_without_nan, :);
%     all_windowed_mean_PS_vectors = data_file.all_windowed_mean_PS_vectors(:,:,:);
%     all_windowed_mean_PS_vectors = all_windowed_mean_PS_vectors(:, :, rows_without_nan);
%     % Reshape the matrix to combine the first two dimensions (891x40) into one
%     reshaped_matrix = reshape(all_windowed_mean_PS_vectors, [], size(all_windowed_mean_PS_vectors, 3));
%     % Check for NaNs in the reshaped matrix
%     nan_vector = any(isnan(reshaped_matrix), 1);
%     clear reshaped_matrix
%     % Convert the result to a column vector
%     nan_vector = nan_vector(:);
%     % Create a logical array indicating rows to keep
%     % kept_rows = rows_without_nan & ~nan_vector;
%     kept_rows = rows_without_nan;
%     kept_rows(rows_without_nan) = ~nan_vector;
% 
%     % Filter the data based on kept rows
%     % all_window_mean_PS_vectors = all_windowed_mean_PS_vectors(:, :, ~nan_vector);
%     label_table = label_table(~nan_vector, :);
% end