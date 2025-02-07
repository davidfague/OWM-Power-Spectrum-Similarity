clear
clc
% fit glme; anova test avg_diff
% load('all_patients_EMS_vs_EFS_table.mat')
load('../AA_Processed Data/allpatients gammamod allregions allitem enc1 correct/all_patients_table.mat')
sorted_combined_table = all_patients_table;
clear all_patients_table

addpath("../subfunctions")
%% filter out <missing> anat, nan/inf ES avgs,
not_missing = ~ismissing(string(sorted_combined_table.anatomical_label));
sorted_combined_table = sorted_combined_table(not_missing,:);

for i=1:size(sorted_combined_table,1)
    sorted_combined_table.avg_diff(i)=mean(sorted_combined_table.EM_EF_diff(i),2);
end


not_nan = ~isnan(sorted_combined_table.avg_diff);
sorted_combined_table=sorted_combined_table(not_nan,:);

not_inf = ~isinf(sorted_combined_table.avg_diff);
sorted_combined_table = sorted_combined_table(not_inf,:);

%% break up anat_patient column
% anat_patient = char(sorted_combined_table.anatomical_label(:));
% anat_patient = cellstr(anat_patient); % Convert to a cell array if it isn't already
% 
% % Split each entry in the cell array into anats and patients
% split_labels = cellfun(@(x) strsplit(x, '_'), anat_patient, 'UniformOutput', false);
% anats = cellfun(@(x) x{1}, split_labels, 'UniformOutput', false);
% patients = cellfun(@(x) x{2}, split_labels, 'UniformOutput', false);
% 
% sorted_combined_table.patient = patients;
% sorted_combined_table.anatomy = anats;

sorted_combined_table.anatomical_label = string(sorted_combined_table.anatomical_label);

sorted_combined_table.patient = categorical(sorted_combined_table.patient_ID);
sorted_combined_table.anatomy = categorical(sorted_combined_table.anatomical_label);

%% Set categories to have 'R Precentral Gyrus' as the baseline for analysis
baseline_category = {'R Precentral Gyrus'};
all_categories = categories(sorted_combined_table.anatomy);
remaining_categories = setdiff(all_categories, baseline_category, 'stable');
new_order = [baseline_category; remaining_categories];
sorted_combined_table.anatomy = reordercats(sorted_combined_table.anatomy, new_order);

clearvars -except sorted_combined_table
%% fit diff lme to the table.
lme = fitlme(sorted_combined_table, ...
    'avg_diff ~ anatomy + (1|channel_ID:patient)', ...
    'FitMethod', 'REML');

lme = fitlme(sorted_combined_table, ...
    'avg_diff ~ anatomy + (1|patient)', ...
    'FitMethod', 'REML');

lme = fitlme(sorted_combined_table, ...
    'avg_diff ~ anatomy + (1|trial_ID:patient)', ...
    'FitMethod', 'REML');

% lme = fitlme(sorted_combined_table, ...
%     'Average_EMS ~ anatomy + (1|patient)', ...
%     'FitMethod', 'REML');


%% check results and save
anova_results = anova(lme);
disp(anova_results);

coeff_results_sorted_by_pValue = sortrows(lme.Coefficients, 'pValue', 'ascend');
disp(coeff_results_sorted_by_pValue)

coeff_results_sorted_by_estimate = sortrows(lme.Coefficients, 'Estimate', 'descend');
disp(coeff_results_sorted_by_estimate)

save("glme_stats.mat", "lme", "anova_results", "coeff_results_sorted_by_estimate", "coeff_results_sorted_by_pValue")

%% correct vs wrong; avg_diff
% correct: 3 correct
% incorrect: > 1 incorrect

sorted_combined_table = compute_correctness_group(sorted_combined_table);
[difference_table] = compute_diff_btwn_corr_incorr(sorted_combined_table);

function [label_table] = compute_correctness_group(label_table)
    label_table.correctness_group = strings(height(label_table), 1);
    
    % Define the new trial correctness classes based on the updated criteria
    for i = 1:height(label_table)
        if sum(label_table.encoding_correctness(i, :), 2) == 3
            label_table.correctness_group(i) = "Correct_Trial";
        elseif sum(label_table.encoding_correctness(i, :), 2) < 2
            label_table.correctness_group(i) = "Incorrect_Trial";
        else
            label_table.correctness_group(i) = "Other";
        end
    end
end

function [difference_table] = compute_diff_btwn_corr_incorr(sorted_table)

% Initialize a new table to store the results
unique_labels = unique(sorted_table.anatomical_label);
n_labels = numel(unique_labels);

% Preallocate the new table
difference_table = table('Size', [n_labels, 4], ...
                         'VariableTypes', {'string', 'double', 'double', 'double'}, ...
                         'VariableNames', {'anatomical_label', 'avg_diff_Correct', 'avg_diff_Incorrect', 'avg_diff_Difference'});

% Loop through each unique anatomical label
for i = 1:n_labels
    % Get the current label
    current_label = unique_labels{i};
    
    % Find the rows corresponding to this label and each correctness group
    correct_row = strcmp(sorted_table.anatomical_label, current_label) & strcmp(sorted_table.correctness_group, "Correct_Trial");
    incorrect_row = strcmp(sorted_table.anatomical_label, current_label) & strcmp(sorted_table.correctness_group, "Incorrect_Trial");
    
    % Extract the avg_diff for correct and incorrect trials
    avg_diff_correct = sorted_table.avg_diff(correct_row);
    avg_diff_incorrect = sorted_table.avg_diff(incorrect_row);
    
    % Store the values in the new table
    difference_table.anatomical_label(i) = current_label;
    difference_table.avg_diff_Correct(i) = avg_diff_correct;
    difference_table.avg_diff_Incorrect(i) = avg_diff_incorrect;
    difference_table.avg_diff_Difference(i) = avg_diff_correct - avg_diff_incorrect;
end

% Calculate the percent change from correct
difference_table.Percent_Change_From_Correct = (difference_table.avg_diff_Difference ./ difference_table.avg_diff_Correct) * 100;

% Display the resulting table
disp(difference_table);
end
