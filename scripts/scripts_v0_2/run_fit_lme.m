clear
clc
% fit glme; anova test avg_diff
load('all_patients_EMS_vs_EFS_table.mat')

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
anat_patient = char(sorted_combined_table.anatomical_label(:));
anat_patient = cellstr(anat_patient); % Convert to a cell array if it isn't already

% Split each entry in the cell array into anats and patients
split_labels = cellfun(@(x) strsplit(x, '_'), anat_patient, 'UniformOutput', false);
anats = cellfun(@(x) x{1}, split_labels, 'UniformOutput', false);
patients = cellfun(@(x) x{2}, split_labels, 'UniformOutput', false);

sorted_combined_table.patient = patients;
sorted_combined_table.anatomy = anats;

sorted_combined_table.patient = categorical(sorted_combined_table.patient);
sorted_combined_table.anatomy = categorical(sorted_combined_table.anatomy);

%% Set categories to have 'R Precentral Gyrus' as the baseline for analysis
baseline_category = {'R Precentral Gyrus'};
all_categories = categories(sorted_combined_table.anatomy);
remaining_categories = setdiff(all_categories, baseline_category, 'stable');
new_order = [baseline_category; remaining_categories];
sorted_combined_table.anatomy = reordercats(sorted_combined_table.anatomy, new_order);

clearvars -except sorted_combined_table
%% fit diff lme to the table.
lme = fitlme(sorted_combined_table, ...
    'avg_diff ~ anatomy + (1|GroupCount:patient)', ...
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