load("summary_p_tables_rebase_all_patients.mat")
all_patients_p_table = sortrows(all_patients_p_table, "p", "descend");
%% remove missing
all_patients_p_table = all_patients_p_table(~ismissing(all_patients_p_table.anat),:);

%% make all lower case 
all_patients_p_table.anat = lower(all_patients_p_table.anat);

%% make column for merge Left and Right
all_patients_p_table.anat_merged = removeDirectionPrefix(all_patients_p_table.anat);

all_patients_p_table = all_patients_p_table(~isnan(all_patients_p_table.p),:);