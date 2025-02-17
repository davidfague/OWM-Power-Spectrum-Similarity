custom_params.hellbender = false;
params = get_parameters(custom_params);

%%
patient_to_images = containers.Map("KeyType", "char", "ValueType", "any");

% gather the PSVs for this region across all patients, channels
for k12wm = [true, false]
    custom_params.k12wm = k12wm;
    params = get_parameters(custom_params);

    for pat_idx = 1:length(params.patient_IDs)
        patient_ID = params.patient_IDs(pat_idx);

        patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, patient_ID);

        for session_idx = 1%:length(patient_preprocessed_data_paths)
            params.session_id = session_idx;
            patient_preprocessed_data_path = patient_preprocessed_data_paths{session_idx};
        
            if params.k12wm
                if params.hellbender
                    prefix = strsplit(patient_preprocessed_data_path,'/');
                else
                    prefix = strsplit(patient_preprocessed_data_path,'\');
                end
                prefix = prefix{end};
                OWM_trial_info_file = load(fullfile(patient_preprocessed_data_path, sprintf("%s_OWM_trialinfo.mat", prefix)));
            
            else
                OWM_trial_info_file = matfile(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"));
            end
        
            images = OWM_trial_info_file.C;
            patient_to_images(sprintf("%d_%d", patient_ID, session_idx)) = images;
        end
    end
end

%% compute likeness % NOTE that the 1 and 2 here are just seperate iterators and not sessions.

% Get all patient/session keys
patients_sessions = patient_to_images.keys;
n_sessions = length(patients_sessions);

% Preallocate a 4D logical array.
% Dimensions: [patient1, patient2, image index in patient1, image index in patient2]
image_likeness = nan(n_sessions, n_sessions, 9, 9);

for idx_1 = 1:n_sessions
    patient_session_1 = patients_sessions{idx_1};
    images_for_patient_session_1 = patient_to_images(patient_session_1);

    for idx_2 = 1:n_sessions
        patient_session_2 = patients_sessions{idx_2};
        images_for_patient_session_2 = patient_to_images(patient_session_2);

        for image_idx1 = 1:9
            image1 = images_for_patient_session_1(image_idx1);
            for image_idx2 = 1:9
                image2 = images_for_patient_session_2(image_idx2);

                if strcmp(image1{1}, image2{1})
                    image_likeness(patient_session_1, patient_session_2, image_idx1, image_idx2) = true;
                else
                    image_likeness(patient_session_1, patient_session_2, image_idx1, image_idx2) = false;
                end
            end
        end
    end
end

%% visualize how many patients share a given image.

% We use a containers.Map to count the number of patients that include a given image.
patient_counts = containers.Map();

for idx = 1:n_sessions
    patient_session = patients_sessions{idx};
    images_for_patient = patient_to_images(patient_session);
    
    % Get the unique images for this patient (to avoid double counting)
    % If your images are stored directly as strings:
    unique_images = unique(images_for_patient);
    % If instead they are stored as cells (like { 'imageName' }), then use:
    % unique_images = unique(cellfun(@(img) img{1}, images_for_patient, 'UniformOutput', false));
    
    for j = 1:length(unique_images)
        imgName = unique_images{j};
        if isKey(patient_counts, imgName)
            patient_counts(imgName) = patient_counts(imgName) + 1;
        else
            patient_counts(imgName) = 1;
        end
    end
end

% Display the count for each unique image.
unique_image_keys = patient_counts.keys;
for i = 1:length(unique_image_keys)
    fprintf('%s is in %d patients.\n', unique_image_keys{i}, patient_counts(unique_image_keys{i}));
end