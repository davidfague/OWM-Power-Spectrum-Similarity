
custom_params.k12wm = false;
custom_params.output_folder_name = 'middle_fixation_baseline';
custom_params.hellbender = false;
params = get_parameters(custom_params);

plot_params = struct();
plot_params.patient_id = 201901;
plot_params.chan_id = 2;
plot_params.image_id = 8;
plot_params.enc_window_ids = params.enc1_win_IDs; % change with enc_id
plot_params.enc_id = 1;
params.session_id = 1;

% for PSVs
plot_params.frequencies_to_use = 1:40; % empty for whatever the data is % for plotting PSVs only

es_freq_bands = {[1:8], [8:20], [20:40]};

%% get all_label_table

all_label_table = table();

% gather label tables
for k12wm = [true, false]
    custom_params.k12wm = k12wm;
    params = get_parameters(custom_params);

    for pat_idx = 1:length(params.patient_IDs)
        patient_id = params.patient_IDs(pat_idx);

        if params.k12wm
            session_ids = get_available_session_ids(params, patient_id);
        else
            session_ids = 1;
        end

        for session_id = session_ids
            params.session_id = session_id;
    
            PS_file = get_PS_file(params, patient_id, false);
            load(PS_file, "label_table")
            all_label_table = [all_label_table; label_table];

        end
    end
end

disp(unique(all_label_table.patient_ID));
disp(unique(all_label_table.session_ID));

save("all_label_table.mat", "all_label_table");
%% gather PSVs for this anat across patients.
% problem: image_id -> image is different across patients.
es_freq_band = [8:30];
% do all images first
possible_PSVs = nan([681, length(es_freq_band), sum(sum(all_label_table.encoding_correctness(:,:), 2)==3)]);

last_row = 1;

% gather the PSVs for this region across all patients, channels
for k12wm = [true, false]
    custom_params.k12wm = k12wm;
    params = get_parameters(custom_params);

    for pat_idx = 1:length(params.patient_IDs)
        patient_id = params.patient_IDs(pat_idx);

        if params.k12wm
            session_ids = get_available_session_ids(params, patient_id);
        else
            session_ids = 1;
        end

        for session_id = session_ids
            params.session_id = session_id;
            fprintf("gathering PSVs for p%d s%d\n", patient_id, session_id)
    
            PS_file = get_PS_file(params, patient_id, false);
            load(PS_file, "label_table")
            rows_to_use = sum(label_table.encoding_correctness(:,:), 2)==3;
            clear label_table
            
            next_last_row = last_row + sum(rows_to_use) - 1;

            load(PS_file, "all_windowed_mean_PS_vectors")

            if size(all_windowed_mean_PS_vectors, 1) == 881
                all_windowed_mean_PS_vectors = all_windowed_mean_PS_vectors(1:681,:,:);
            end

            if length(last_row:next_last_row) ~= sum(rows_to_use)
                error("length(last_row:next_last_row) ~= sum(rows_to_use): %d %d", length(last_row:next_last_row), sum(rows_to_use))
            elseif length(es_freq_band) ~= size(possible_PSVs, 2)
                error("length(es_freq_band) ~= size(possible_PSVs, 2): %d %d", length(es_freq_band), size(possible_PSVs, 2))
            elseif size(all_windowed_mean_PS_vectors,1) ~= size(possible_PSVs, 1)
                error("size(all_windowed_mean_PS_vectors,1) ~= size(possible_PSVs, 1): %d %d", size(all_windowed_mean_PS_vectors,1), size(possible_PSVs, 1))
            end

            possible_PSVs(:,:,last_row:next_last_row) = all_windowed_mean_PS_vectors(:, es_freq_band, rows_to_use);
            clear all_windowed_mean_PS_vectors
            last_row = next_last_row;
        end
    end
end
plot_params.frequencies_to_use = [];
plot_params.anat = "All-Anat";
title_str = sprintf("Correct-Trials All-Patients %d-%d", min(es_freq_band), max(es_freq_band));
plot_params.enc_window_ids = params.enc1_win_IDs;
enforced_clim = [];
fig = plot_PSV(mean(possible_PSVs(:,params.maint_win_IDs,:), [2 3], 'omitnan'), plot_params, title_str, enforced_clim);

%% this is channel-specific method (not used here)

% plot_PSVs(params, plot_params, false) % plot trials with the image
% close all
% 
% % BI
% plot_PSVs(params, plot_params, true) % plot trials without the image
% close all