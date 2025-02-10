%% parameters
analysis_params = struct();
analysis_params.min_npatient_per_region = 3;

p_threshhold = 0.975;
%% load data
% load('summary_utah_patients_table.mat')
all = load("../results/new results/summary_p_tables_rebase_all_patients.mat", "all_patients_p_table");
utah = load("../results/new results/summary_p_tables_rebase_utah_all.mat", "all_patients_p_table");
k12wm = load('../results/new results/summary_p_tables_rebase_k12wm_all.mat', "all_patients_p_table");

all_patients_p_table = all.all_patients_p_table; % select table

% addpath("../subfunctions/")
params = get_parameters();

%% remove missing
all_patients_p_table = all_patients_p_table(~ismissing(all_patients_p_table.anat),:);

%% make all lower case 
all_patients_p_table.anat = lower(all_patients_p_table.anat);

%% make column for merge Left and Right
all_patients_p_table.anat_merged = removeDirectionPrefix(all_patients_p_table.anat);

%% generate new table with total occurences of the anat, percent of observations that are > p_threshold
% Get the list of unique anatomical regions (using the merged names)
regions = unique(all_patients_p_table.anat_merged);

% Preallocate arrays to store summary statistics for each region
nRegions = numel(regions);
totalOccurrences = zeros(nRegions, 1);
fracAboveThreshold = zeros(nRegions, 1);
totalAboveThreshold = zeros(nRegions, 1);
nUniquePatients = zeros(nRegions, 1);

% Loop over each region to compute the stats
for i = 1:nRegions
    % Identify rows corresponding to the current region
    idx = strcmp(all_patients_p_table.anat_merged, regions{i});
    
    % Total number of observations in this region
    totalOccurrences(i) = sum(idx);
    
    % Fraction of observations with p > p_threshhold
    fracAboveThreshold(i) = sum(all_patients_p_table.p(idx) > p_threshhold) / totalOccurrences(i);
    totalAboveThreshold(i) = sum(all_patients_p_table.p(idx) > p_threshhold);
    
    % Count unique patients for this region (assuming the column is named 'patient_id')
    nUniquePatients(i) = numel(unique(all_patients_p_table.patient_id(idx)));
end

% Create the summary table
summaryTable = table(regions, totalOccurrences, totalAboveThreshold, fracAboveThreshold, nUniquePatients, ...
    'VariableNames', {'anat_merged', 'total_occurrences', 'total_above_threshold', 'frac_above_threshold', 'n_unique_patients'});

% % Filter out regions that do not have at least the minimum number of unique patients
% summaryTable = summaryTable(summaryTable.n_unique_patients >= analysis_params.min_npatient_per_region, :);

% % Display the summary table
% disp(summaryTable);

%% filter region

load("../results/new results/regions_to_use.mat", "regions_to_use")
summaryTable_unfiltered = summaryTable;

% filter based on regions meeting criteria in the all_patient_table
summaryTable = summaryTable(ismember(summaryTable.anat_merged, regions_to_use), :);

% add placeholder for missing regions
for region_idx = 1:length(regions_to_use)
    region = regions_to_use(region_idx);
    if ~ismember(region, summaryTable.anat_merged)
        summaryTable(end,:) = table(region, 0, nan, nan, 0, ...
            'VariableNames', {'anat_merged', 'total_occurrences', 'total_above_threshold', 'frac_above_threshold', 'n_unique_patients'});
    end
end

%%

% Example MATLAB code to create the desired bar graph

% --- Assume your data is already loaded into the table:
% summary_p_tables_work_all_patients
% with fields: anat_merged, total_occurrences, frac_above_threshold, n_unique_patients

% For example, if loading from a CSV file:
% summary_p_tables_work_all_patients = readtable('your_data.csv');

% Sort the table by descending frac_above_threshold
[sortedFrac, sortIdx] = sort(summaryTable.frac_above_threshold, 'descend');
sortedAnat   = summaryTable.anat_merged(sortIdx);
sortedOcc    = summaryTable.total_occurrences(sortIdx);
sortedTotal    = summaryTable.total_above_threshold(sortIdx);

% Create new labels that include the anatomical region and the occurrence count (N)
labels = cell(length(sortedAnat), 1);
for i = 1:length(sortedAnat)
    labels{i} = sprintf('%s (N=%d)', sortedAnat{i}, sortedOcc(i));
end

% Create a horizontal bar graph
figure('WindowState','maximized');
barh(sortedFrac, 'FaceColor', [0.2, 0.6, 0.8]);  % You can adjust the color as desired

% By default, barh places the first element at the bottom. To display the highest value on top, reverse the y-axis:
set(gca, 'YDir', 'reverse');

% Set the y-axis tick labels to our custom labels
set(gca, 'YTick', 1:length(labels), 'YTickLabel', labels);

% Add axis labels and a title
xlabel('Fraction Above Threshold');
ylabel('Anatomical Region (with Occurrences)');
title('Anatomical Regions Sorted by Descending Fraction Above Threshold');
% xlim([0 0.5]);

% Optional: Adjust font size for better readability
set(gca, 'FontSize', 12);
