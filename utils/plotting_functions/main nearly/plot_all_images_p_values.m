function plot_all_images_p_values(params)
% This function gathers p-values and anatomical labels for channels across patients
% and plots sorted p-values with corresponding channel anatomical labels.

    % Initialization
    enc_id = 1; % Encoding ID (can be adjusted based on requirements)
    all_p = []; % To store concatenated p-values
    all_anat = {}; % To store concatenated anatomical labels (as cell array)

    % Gather data across patients
    for patient_id = params.patient_IDs
        % Load data for the current patient
        loaded_data = load(sprintf("results/t_testing/%s_%d.mat", patient_id, enc_id));
        
        % Concatenate p-values and anatomical labels
        all_p = [all_p; loaded_data.p_values_by_chan]; % Append rows of p-values
        all_anat = [all_anat; loaded_data.anat_by_chan]; % Append anatomical labels
    end

    % Ensure all_anat is a column cell array
    all_anat = all_anat(:);

    % Sort data by descending order of the first p-value column (or choose a specific criterion)
    [sorted_p, sort_idx] = sort(all_p(:, 1), 'descend'); % Sorting based on the first column of p-values
    sorted_anat = all_anat(sort_idx);

    % Plot sorted data
    figure;
    bar(sorted_p)
    xticks(1:length(sorted_anat)); % Ensure ticks correspond to the sorted data
    xticklabels(sorted_anat) % Use anatomical labels as tick labels
    ylabel('p-value')
    xlabel('Channels')
    title('Sorted p-values by Channel')
    
    % Rotate x-axis labels for better visibility
    xtickangle(90)

    % Add a horizontal reference line (e.g., significance threshold)
    hold on;
    yline(0.05, '--r', 'Significance Threshold');
    hold off;
end
% 
% function plot_all_images_p_values(params)
% % WIP code for plotting after finishing computing.
%     enc_id = 1;
%     all_p = [];
%     all_anat = [];
%     for patient_id = params.patient_IDs % gather data across patients
% 
%         loaded_data = load(sprintf("results/t_testing/%s_%s.mat", patient_id, enc_id));
%         all_p = [all_p, loaded_data.p_values_by_chan]; % loaded_data.p_values_by_chan is size (n_channel_ids_to_use, n_image_ids=9) # n_channel_ids_to_use can vary by patient;
%         all_anat = [all_anat, loaded_data.anat_by_chan]; % loaded_data.p_values_by_chan is size (n_channel_ids_to_use)
%         % loaded_data.chan_id_by_chan
%     end
% 
%     % make all_p to be size(sum_all_n_channel_ids_to_use, n_image_ids=9)
%     % make all_anat to be size(sum_all_n_channel_ids_to_use)
% 
%     % plot
% 
%     figure;
%     % Sort by descending order of final_p_by_channels
%     [sorted_p, sort_idx] = sort(all_p, 'descend');
%     sorted_anat = all_anat(sort_idx);
% 
%     % Plot the sorted data
%     bar(sorted_p)
%     xticks(1:length(sorted_anat)); % Ensure all ticks are present
%     xticklabels(sorted_anat)
%     ylabel('p-value')
%     xlabel('Channels')
%     title('Sorted p-values by Channel')
% 
%     % Rotate x-axis labels for better visibility
%     xtickangle(90)
% 
% end