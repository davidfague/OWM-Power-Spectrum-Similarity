clear
clc
% fit glme; anova test avg_diff
% load('all_patients_EMS_vs_EFS_table.mat')
load('all_patients_table.mat')
sorted_combined_table=all_patients_table;

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


% Find and replace in the '' column
sorted_combined_table.anatomy(strcmp(sorted_combined_table.anatomy, findString)) = {replaceString};
    
% remove certain labels

sorted_combined_table(sorted_combined_table.anatomy == 'WM', :) = [];
sorted_combined_table(sorted_combined_table.anatomy == 'Ventricle', :) = [];
sorted_combined_table(sorted_combined_table.anatomy == 'Unknown', :) = [];
sorted_combined_table(sorted_combined_table.anatomy == 'Resection cavity', :) = [];

sorted_combined_table(sorted_combined_table.anatomy == 'Left Cerebral White Matter', :) = [];

sorted_combined_table(sorted_combined_table.anatomy == 'Right Cerebral White Matter', :) = [];

tabulatedrois = tabulate(sorted_combined_table.anatomy);
lowrois=tabulatedrois;
lowrois=cell2table(lowrois);
lowrois_labels=lowrois.lowrois1(lowrois.lowrois2<80);

for i=1:length(lowrois_labels)
    sorted_combined_table(sorted_combined_table.anatomy == lowrois_labels(i), :) = [];
end
sorted_combined_table.anatomy = removecats(sorted_combined_table.anatomy, lowrois_labels);

[C,ia,ic] = unique(sorted_combined_table.anatomy)
tabulatedrois = tabulate(sorted_combined_table.anatomy);



%% Set categories to have 'R Precentral Gyrus' as the baseline for analysis
baseline_category = {'L Precentral Gyrus'};
all_categories = categories(sorted_combined_table.anatomy);
remaining_categories = setdiff(all_categories, baseline_category, 'stable');
new_order = [baseline_category; remaining_categories];
sorted_combined_table.anatomy = reordercats(sorted_combined_table.anatomy, new_order);
figure; histogram(sorted_combined_table.anatomy)
clearvars -except sorted_combined_table
%% fit diff lme to the table.
lme = fitlme(sorted_combined_table, ...
    'avg_diff ~ anatomy + (1|channel_ID:patient)', ...
    'FitMethod', 'REML');

lme = fitlme(sorted_combined_table, ...
    'avg_diff ~ anatomy + (1|patient)', ...
    'FitMethod', 'REML');
%%
lme = fitlme(sorted_combined_table, ...
    'avg_diff ~ anatomy + (1|trial_ID:patient)', ...
    'FitMethod', 'REML');


% check results and save
anova_results = anova(lme);
disp(anova_results);

coeff_results_sorted_by_pValue = sortrows(lme.Coefficients, 'pValue', 'ascend');
disp(coeff_results_sorted_by_pValue)

clear siglabels
siglabels=table();

for i=2:18
     chr=string(coeff_results_sorted_by_pValue{i,1});
    newChr = strrep(chr,'anatomy_','')
siglabels.lab(i-1)=newChr;
end

clear rmsiglabels
rmsiglabels=table();
count=1;
for i=19:52
    chr=string(coeff_results_sorted_by_pValue{i,1});
    newChr = strrep(chr,'anatomy_','')
rmsiglabels.lab(count)=newChr;
count=count+1;
end


for i=1:height(rmsiglabels)

    sorted_combined_table(sorted_combined_table.anatomy == rmsiglabels.lab(i), :) = [];
end
sorted_combined_table.anatomy = removecats(sorted_combined_table.anatomy, rmsiglabels.lab);

%%
coeff_results_sorted_by_estimate = sortrows(lme.Coefficients, 'Estimate', 'descend');
disp(coeff_results_sorted_by_estimate)
%%
save("glme_stats.mat", "lme", "anova_results", "coeff_results_sorted_by_estimate", "coeff_results_sorted_by_pValue")

%% plot



%% plot specific ROIs


