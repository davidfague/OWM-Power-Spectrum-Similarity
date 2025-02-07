all = load("summary_p_tables_rebase_all_patients.mat", "all_patients_p_table");
utah = load("summary_p_tables_rebase_utah_all.mat", "all_patients_p_table");
k12wm = load('summary_p_tables_rebase_k12wm_all.mat', "all_patients_p_table");

ismember(unique(utah.all_patients_p_table.anat), unique(all.all_patients_p_table.anat))