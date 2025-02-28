% scripts/
%% get parameters
script_specs = struct();
script_specs.clip_inf_similarities = true;
script_specs.skip_existing = false;

custom_params = get_custom_params(script_specs); % get other defaults like hellbender, ouput folder

params = get_parameters(custom_params);

clear script_specs custom_params

%%
params.patient_id = 201910;
params.chan_id = 66;
params.image_id = 6;
params.enc_window_ids = params.enc1_win_IDs; % change with enc_id
params.enc_id = 1;
params.session_id = 1;

params.freq_min = 1;
params.freq_max = 40;

% for PSVs
params.frequencies_to_use = 1:40; % for plotting PSVs only % empty for whatever the data is 

% for WI vs BI
params.type = 'EMS'; % for plotting WI_vs_BI only
params.mean_out_time_dimensions = true;
params.average_diff = true;
params.same_n = true; % only affects plot_similarity_means_heatmaps

params.only_all3_correct = true; % filters test trials (with item)

% es_freq_bands = {[1:8], [8:20], [20:40], [1:40]};
params.es_freq_bands = {[1:40]};%{[4:8],[8:12],[12:30], [30:70]};
% es_freq_bands = {[1:40]};
% look at power spectra to determine bands
% 4:8
% 8:12-14
% 15:30
% 30:70
% 70:140

%% params that will need updated if their input is changed below
% might need to put the above params in script_specs and then put this
% logic toward end of get_parameters and use 'if the required fields are present'.

params.patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, params.patient_id);
params.patient_preprocessed_data_path = params.patient_preprocessed_data_paths{params.session_id};

params.WI_BI_folder_to_save_in = fullfile(sprintf("results/WI vs BI/p%s chan%s image%s enc%s sess%d", ...
            num2str(params.patient_id), num2str(params.chan_id), ...
            num2str(params.image_id), num2str(params.enc_id), params.session_id));

params.anat_labels = get_anat_labels(params.patient_preprocessed_data_path, params);

% params.image_labels = load(fullfile(params.patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted

params.anat = string(params.anat_labels.labelsanatbkedit.anatmacro1(params.chan_id));
%% compute table if needed (can probably have a seperate script for this)
params.recompute_table = true;
% params.patient_IDs = [010];

if params.recompute_table
    original_mean_out_time = params.mean_out_time_dimensions;
    params.mean_out_time_dimensions = true;
    for freqs_idx = 1:length(params.es_freq_bands)
        all_p_table = table();
        params.ES_freq_band = params.es_freq_bands{freqs_idx};
        params.freq_min = min(params.ES_freq_band);
        params.freq_max = max(params.ES_freq_band);

        for pat_idx = 1:length(params.patient_IDs)
            params.patient_id = params.patient_IDs(pat_idx);
            all_p_table = calc_WI_vs_BI_table(params, all_p_table);
        end

        if params.k12wm
            suffix = 'k12wm';
        else
            suffix = 'utah';
        end

        save_file = sprintf('summary_p_tables_midbase_%d-%dhz_%s_clip_infs.mat', params.freq_min, params.freq_max, suffix);
        fprintf("saving %s", save_file)

        save(save_file, 'all_p_table')
    end
end
params.mean_out_time_dimensions = original_mean_out_time;
clearvars patient_ids_to_use table_all_patients original_mean_out_time

%% compute significance and clustering. (plot_WI_vs_BI does the significance calculating and plotting, calling some functions that calc_WI_vs_BI_table also calls. Probably need to modularize the significance calculating functions.)
% % this has a separate goal from the table, but can be used to iterate
% % across the table. some other script was used to get "data" subset from a
% % table comuted above for example.

% params.patient_id = 201908;
% params.image_id = 7;
% params.chan_id = 17;
% for row_idx = 1:length(data.patient_id)
%     params.chan_id = data(row_idx,:).chan_id;
%     params.image_id = data(row_idx,:).image_id;
%     params.patient_id = data(row_idx,:).patient_id;
    % fprintf("row%s p%s chan%s image%s\n",num2str(row_idx), num2str(params.patient_id), num2str(params.chan_id), num2str(params.image_id))
    % % plot PSVs
    % WI
    plot_PSVs(params, false); % plot trials with the image
    % close all
    
    % BI
    plot_PSVs(params, true); % plot trials without the image
    % close all
    
    % % plot WI vs BI EMS, t-test, clusters
    % plot WI vs BI EMSs
    % for ave_diff = [true, false]
    
        % params.average_diff = ave_diff;
    WI_vs_BI(params, true);
    close all
% end

%% plot within_trial
% params.patient_id = 201908;
% params.image_id = 7;
% params.chan_id = 17;
% params.without_image = false;
% plot_within_trial(params, params)
% close all

% for row_idx = 8:length(data.patient_id)
    % params.chan_id = data(row_idx,:).chan_id;
    % params.image_id = data(row_idx,:).image_id;
    % params.patient_id = data(row_idx,:).patient_id;

    % with image
    for without_image = [true, false]
        params.without_image = without_image;
        plot_within_trial(params)
        close all
    end

% end

%% functions
% have been stored away in utils.

%% legacy (now stored elsewhere)
% function plot_PSVs(params, params, without_image)
%     folder_to_save_in = sprintf("results/PSV plots/same_clims/p%s chan%s image%s enc%s", ...
%                 num2str(params.patient_id), num2str(params.chan_id), ...
%                 num2str(params.image_id), num2str(params.enc_id));
%     mkdir(folder_to_save_in)
%     mkdir(sprintf("%s/single obs/", folder_to_save_in))
%     mkdir(sprintf("%s/mean across trials/", folder_to_save_in))
% 
%     % load PSV
%     PS_file = get_PS_file(params.output_folder, params.patient_id, false);
%     load(PS_file, "label_table")
% 
%     % subset table by params
%     % rows with the image and params
%     if ~without_image
%         rows_to_use = label_table.channel_ID == params.chan_id & ... % channel
%         label_table.patient_ID == params.patient_id & ... % patient
%         label_table.encID_to_imageID(:,params.enc_id)==params.image_id & ... % image in first encoding
%         sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
%     else
%         % rows without the image
%         rows_to_use = label_table.channel_ID == params.chan_id & ... % channel
%         label_table.patient_ID == params.patient_id & ... % patient
%         sum(label_table.encID_to_imageID(:,:)~=params.image_id,2)==3 & ... % image in first encoding
%         sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
%     end
%     subset_table = label_table(rows_to_use,:); % can check
% 
%     % get anat from table
%     params.anat = unique(string(label_table(rows_to_use,:).anatomical_label)); % get anat from table
%     % clear label_table %now would be the time to do it but probably keep for
%     % using trial_id
%     if length(params.anat) > 1
%         error("should only detect one anat")
%     end
% 
%     % get correspoding rows (trials) of PSVs from PS_file
%     % PS_file = matfile(PS_file);
%     load(PS_file, 'all_windowed_mean_PS_vectors')
%     possible_PSVs_for_this_obs = all_windowed_mean_PS_vectors(:,:,rows_to_use);
%     rows_without_nans =~any(isnan(possible_PSVs_for_this_obs),[1 2]);
%     possible_PSVs_for_this_obs = possible_PSVs_for_this_obs(:,:,rows_without_nans); % remove trials with nans
%     subset_table = subset_table(rows_without_nans,:);
%     clear all_windowed_mean_PS_vectors
% 
%     % flattened_data = possible_PSVs_for_this_obs(params.enc_window_ids, :, :);
%     % Calculate the mean and standard deviation of the flattened data
%     % data_mean = mean(flattened_data(:));
%     % data_std = std(flattened_data(:));
%     % Set the color limits to be 2 standard deviations around the mean
%     % clims_to_use = [data_mean - 2*data_std, data_mean + 2*data_std];
% 
%     clims_to_use = [-1, 2];% with:[-3.3583    4.5770]; without: [-3.5836    4.3866]
%     max_iter = 9;
%     max_iter_to_use = min(max_iter, size(possible_PSVs_for_this_obs,3));
%     for obs_slice = 1:max_iter_to_use % can also select 1
%         if ~without_image
%             title_str = sprintf("p%s chan%s image%s enc%s trial%s", ...
%                 num2str(params.patient_id), num2str(params.chan_id), ...
%                 num2str(params.image_id), num2str(params.enc_id), ...
%                 num2str(subset_table(obs_slice,:).trial_ID));
%         else
%              title_str = sprintf("p%s chan%s WITHOUT image%s enc%s trial%s", ...
%                 num2str(params.patient_id), num2str(params.chan_id), ...
%                 num2str(params.image_id), num2str(params.enc_id), ...
%                 num2str(subset_table(obs_slice,:).trial_ID));
%         end
%         fig = plot_PSV(possible_PSVs_for_this_obs(:,:,obs_slice), params, title_str, clims_to_use); % single observation
%         saveas(fig, sprintf("%s/single obs/%s obs%s.png",folder_to_save_in, title_str, num2str(obs_slice)))
%     end
% 
% 
%     % flattened_data = mean(possible_PSVs_for_this_obs(params.enc_window_ids,:,:), 3, 'omitnan');
%     % data_mean = mean(flattened_data(:));
%     % data_std = std(flattened_data(:));
%     % Set the color limits to be 2 standard deviations around the mean
%     % clims_to_use = [data_mean - 2*data_std, data_mean + 2*data_std];
%    % with: [-0.8951    2.1138]; without: [-0.5547    1.3577]
%    clims_to_use = [-1, 2];
% 
%     if ~without_image
%         title_str = sprintf("p%s chan%s image%s enc%s mean across trials", ...
%             num2str(params.patient_id), num2str(params.chan_id), ...
%             num2str(params.image_id), num2str(params.enc_id));
%     else
%         title_str = sprintf("p%s chan%s WITHOUT image%s enc%s mean across trials", ...
%             num2str(params.patient_id), num2str(params.chan_id), ...
%             num2str(params.image_id), num2str(params.enc_id));
%     end
%     mean(possible_PSVs_for_this_obs, 3, 'omitnan')
%     fig = plot_PSV(mean(possible_PSVs_for_this_obs, 3, 'omitnan'), params, title_str, clims_to_use); % mean of all obs
%     saveas(fig, sprintf("%s/mean across trials/%s.png", folder_to_save_in, title_str))
% 
% end