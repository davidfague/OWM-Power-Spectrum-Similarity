%% plotting avg encoding similarity over time
addpath('../plotting functions')

output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
mf = matfile(fullfile(output_folder,'201907all3_ES.mat'));
mf_table = matfile(fullfile(output_folder,'201907.mat'));

% label_table = mf_table.label_table;
% label_table.anatomical_label = string(label_table.anatomical_label);
% label_table = compute_correctness_group(label_table);

% pick better slice. plot for regions.

% size(mf.all3_ES_matrix) % 12180          41         640           3
% (channel*trial, encodingWindows, allWindows, encID)            

enc_id = 3;
% iteration = 1;
% ES_matrix = squeeze(mf.all3_ES_matrix(iteration, :,:,enc_id)); % one
% iteration

% combine all iterations (channels x trials)
rows_without_nans = ~any(any(any(isnan(mf.all3_ES_matrix), 2), 3), 4);
indices_without_nans = find(rows_without_nans);
max_id = max(indices_without_nans);
min_id = min(indices_without_nans);
ES_matrix = mf.all3_ES_matrix(min_id:max_id, :,:,enc_id);
ES_matrix = ES_matrix(indices_without_nans - (min_id-1),:,:,:);
% ES_matrix = squeeze(mean(ES_matrix,1)); % do both dimensions at same time
% so std is correct
%% TODO:
% subset mf_table with the same indices
label_table = mf_table.label_table;
label_table = label_table(rows_without_nans,:);
label_table.anatomical_label = string(label_table.anatomical_label);
% subset table and ES_matrix by new anatomy, performance slice.
target_anat = 'L Hippocampus';
plot_avg_ES_over_time_anat(ES_matrix, label_table, target_anat)

% subset by encoding correct
function [ES_matrix, label_table] = subset_by_enc(ES_matrix, label_table, target_enc_ID, target_enc_correctness)
rows_with_enc_target = label_table.encoding_correctness(:, target_enc_ID) == target_enc_correctness;
ES_matrix = ES_matrix(rows_with_enc_target,:,:,:);
label_table = label_table(rows_with_enc_target,:);
end

%%
% clearvars -except ES_matrix
% plot_avg_ES_over_time(ES_matrix)

function plot_avg_ES_over_time_anat(ES_matrix, label_table, target_anat)
rows_with_anat = label_table.anatomical_label == target_anat;
anat_table = label_table(rows_with_anat,:);
anat_matrix = ES_matrix(rows_with_anat,:,:,:);
[anat_matrix, anat_table] = subset_by_enc(anat_matrix, anat_table, 3, 0);
plot_avg_ES_over_time(anat_matrix)
title(['Encoding similarity during all trials for ', target_anat]);
end

function plot_avg_ES_over_time(ES_matrix)
avg_ES_over_time = squeeze(mean(ES_matrix, [1 2]));
% avg_ES_over_time(isinf(avg_ES_over_time)) = 0;
ntime = 1:length(avg_ES_over_time); % can update to do the mean of window start and end time for this windowID.

std_ES_over_time = squeeze(std(ES_matrix, 0, [1 2]));
clear ES_matrix
std_ES_over_time(isnan(std_ES_over_time)) = 0;

figure()
plot(ntime, avg_ES_over_time, 'LineWidth', 2); % Plot the average line
y_limits = ylim;
avg_ES_over_time(isinf(avg_ES_over_time)) = y_limits(2);

upper_bound = avg_ES_over_time + std_ES_over_time;
lower_bound = avg_ES_over_time - std_ES_over_time;

plot(ntime, avg_ES_over_time, 'LineWidth', 2); % Plot the average line
hold on;

% Plot the shaded area using the fill function
fill([ntime, fliplr(ntime)], [upper_bound', fliplr(lower_bound')], ...
     'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Customize plot
xlabel('Time Window'); % change to time
ylabel('Encoding Similarity');
title('Encoding similarity during all trials');
grid on;
add_event_lines()
plot_horizontal_means(avg_ES_over_time)
hold off;
end

% avg_ES_over_time = squeeze(mean(ES_matrix, [1 2]));
% % avg_ES_over_time(isinf(avg_ES_over_time)) = 0;
% ntime = 1:length(avg_ES_over_time); % can update to do the mean of window start and end time for this windowID.
% 
% std_ES_over_time = squeeze(std(ES_matrix, 0, [1 2]));
% clear ES_matrix
% std_ES_over_time(isnan(std_ES_over_time)) = 0;
% 
% figure()
% plot(ntime, avg_ES_over_time, 'LineWidth', 2); % Plot the average line
% y_limits = ylim;
% avg_ES_over_time(isinf(avg_ES_over_time)) = y_limits(2);
% 
% upper_bound = avg_ES_over_time + std_ES_over_time;
% lower_bound = avg_ES_over_time - std_ES_over_time;
% 
% plot(ntime, avg_ES_over_time, 'LineWidth', 2); % Plot the average line
% hold on;
% 
% % Plot the shaded area using the fill function
% fill([ntime, fliplr(ntime)], [upper_bound', fliplr(lower_bound')], ...
%      'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
% % Customize plot
% xlabel('Time Window'); % change to time
% ylabel('Encoding Similarity');
% title('Encoding similarity during the trial');
% grid on;
% add_event_lines()
% plot_horizontal_means(avg_ES_over_time)
% hold off;