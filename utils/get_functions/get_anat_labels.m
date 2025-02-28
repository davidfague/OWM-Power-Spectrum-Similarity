function anat_labels = get_anat_labels(patient_preprocessed_data_path, params)
    if params.k12wm
        if params.hellbender
            prefix = strsplit(patient_preprocessed_data_path,'/'); % linux
        else
            prefix = strsplit(patient_preprocessed_data_path,'\');
        end
        prefix = prefix{end};
        labels = load(fullfile(patient_preprocessed_data_path, sprintf("%s_labelsAnat.mat", prefix)));
        anat_labels = struct();
        anat_labels.labelsanatbkedit = labels.bipolarAnat;
        clear labels prefix
    else
        % image_labels = load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted
        anat_labels = load(fullfile(patient_preprocessed_data_path, ...
        "D_OWM_t_bipolar.mat"), 'labelsanatbkedit'); % i've already made sure this is the same as the table
    end

end