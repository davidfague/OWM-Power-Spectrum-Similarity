load("summary_p_tables_rebase_all_patients.mat", "all_patients_p_table")

addpath("../subfunctions/")

%% remove missing
all_patients_p_table = all_patients_p_table(~ismissing(all_patients_p_table.anat),:);


%% make all lower case 
all_patients_p_table.anat = lower(all_patients_p_table.anat);

%% investigate nan

mean(isnan(all_patients_p_table.p)) % 7.44% nans

all_patients_p_table(isnan(all_patients_p_table.p),:)


%%
regions = unique(all_patients_p_table.anat);

% Preallocate arrays to store summary statistics for each region
nRegions = numel(regions);

% Loop over each region to compute the stats
for i = 1:nRegions
    % Identify rows corresponding to the current region
    idx = strcmp(all_patients_p_table.anat, regions{i});
    
    % Total number of observations in this region
    totalOccurrences(i) = sum(idx);

end