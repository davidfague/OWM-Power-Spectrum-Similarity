%% 
addpath("../subfunctions")
addpath("../plotting functions")
addpath("../plotting functions/main")

params = get_parameters();

%% select data
     % "201915"       42           0    "L Middle Frontal Gyrus"       1 Close p<0.11
     % "201915"       14           0    "L Middle Frontal Gyrus"       2    
     % "201915"       42           0    "L Middle Frontal Gyrus"       3    
     % "201915"       14           0    "L Middle Frontal Gyrus"       4    
     % "201915"       14           0    "L Middle Frontal Gyrus"       5    
     % "201915"       43           0    "L Middle Frontal Gyrus"       5    
     % "201915"       16           0    "L Middle Frontal Gyrus"       6    
     % "201915"       42           0    "L Middle Frontal Gyrus"       6    
     % "201915"       13           0    "L Middle Frontal Gyrus"       7    
     % "201915"       15           0    "L Middle Frontal Gyrus"       7    
     % "201915"       16           0    "L Middle Frontal Gyrus"       7    
     % "201915"       14           0    "L Middle Frontal Gyrus"       8    
     % "201915"       16           0    "L Middle Frontal Gyrus"       8    
     % "201915"       43           0    "L Middle Frontal Gyrus"       8    
     % "201915"       13           0    "L Middle Frontal Gyrus"       9    
     % "201915"       14           0    "L Middle Frontal Gyrus"       9    
     % "201915"       15           0    "L Middle Frontal Gyrus"       9   LESS THAN 0 while null is 0 
plot_params = struct();
plot_params.patient_id = 201915;
plot_params.chan_id = 42;
plot_params.image_id = 3;
plot_params.enc_window_ids = params.enc1_win_IDs; % change with enc_id
plot_params.enc_id = 1;

% for PSVs
plot_params.frequencies_to_use = 1:40; % empty for whatever the data is % for plotting PSVs only

% for WI vs BI
plot_params.type = 'EMS'; % for plotting WI_vs_BI only
plot_params.mean_out_time_dimensions = false;
plot_params.average_diff = true;
plot_params.same_n = true; % only affects plot_similarity_means_heatmaps

%% compute
% plot the PSV
% plot_PSVs(params, plot_params, false) % plot trials with the image
% 
% plot_PSVs(params, plot_params, true) % plot trials without the image
% close all

% plot WI vs BI EMSs
[real_t_values, real_p_values, surrogate_t_values, surrogate_p_values] = plot_WI_vs_BI(params, plot_params);


%% functions
function [real_t_values, real_p_values, surrogate_t_values, surrogate_p_values] = plot_WI_vs_BI(params, plot_params)
    plot_params.WI_BI_folder_to_save_in = sprintf("results/WI vs BI/p%s chan%s image%s enc%s", ...
                num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
                num2str(plot_params.image_id), num2str(plot_params.enc_id));

    folder_to_save_in = plot_params.WI_BI_folder_to_save_in;

    mkdir(folder_to_save_in);

    patient_id = [plot_params.patient_id];
    channel_ids_to_use = [plot_params.chan_id];

    patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
    % image_labels = load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted
    anat_labels = load(fullfile(patient_preprocessed_data_path, ...
    "D_OWM_t_bipolar.mat"), 'labelsanatbkedit'); % i've already made sure this is the same as the table
    
    for chan_idx = 1:length(channel_ids_to_use)
        chan_id = channel_ids_to_use(chan_idx);
        plot_params.anat = string(anat_labels.labelsanatbkedit.anatmacro1(chan_id));
            
        % get WI and BI matrices of sizes: (E_times, M_times, trial_combinations)
        [WI, BI] = get_and_pre_process_WI_BI(plot_params.patient_id, params.comp_options, plot_params.enc_id, ...
            plot_params.image_id, plot_params.chan_id, params, plot_params.type);

        % load previous calculations (not really used)
        % selected_p_data = load(sprintf("results/selected_p_data_p%s_image%s_chan%s.mat", ...
        %     num2str(plot_params.patient_id),num2str(plot_params.image_id),num2str(plot_params.chan_id)));

        % plot EMS traces
        fig = plot_similarity_means_timecourse(WI,BI, plot_params);
        saveas(fig, sprintf("%s/WI_BI_time_courses_same_n.png", folder_to_save_in))

        % plot EMS heatmaps
        for same_n = [true, false]
            plot_params.same_n = same_n;
            fig = plot_similarity_means_heatmap(WI,BI, plot_params);
            if plot_params.same_n
                saveas(fig, sprintf("%s/WI_BI_heatmaps_same_n.png", folder_to_save_in))
            else
                saveas(fig, sprintf("%s/WI_BI_heatmaps.png", folder_to_save_in))
            end
        end

        % temporal generalization
        if plot_params.mean_out_time_dimensions
            WI=mean(WI, [1 2]);
            BI=mean(BI, [1 2]);
        end

        % calculate real diff
        if plot_params.average_diff
            diff = calc_diff_avg(WI, BI, 1000);
        else
            diff = calc_diff_all_pairs(WI, BI, 10000);
        end

        % calculate t-tests on real diff
        [real_t_values, real_p_values] = calc_ttest(diff);
        
        % plot_diff(diff)

        % calculate surrogate and t_test it
        [surrogate_t_values, surrogate_p_values, surrogate_diff] = calc_surrogate(WI, BI, 1000, plot_params.average_diff, plot_params.mean_out_time_dimensions);

        % plot differences
        if plot_params.mean_out_time_dimensions % surrogate_diff is too large without meaning out time dimensions beforehand.
            fig = plot_diff(diff, surrogate_diff(:));
            saveas(fig, sprintf("%s/WI_minus_BI_distributions.png", folder_to_save_in))
        end

        % % plot t_values
        % fig = plot_t_values(real_value, null_values);

        p = cluster_analysis_and_visualization_abs_main(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values, plot_params);
        % final_p = mean(real_p_values > surrogate_p_values);
        final_p = mean(real_t_values > surrogate_t_values);

        % p_by_channels(chan_idx) = final_p;
    end
end


%% legacy (now stored elsewhere)
% function plot_PSVs(params, plot_params, without_image)
%     folder_to_save_in = sprintf("results/PSV plots/same_clims/p%s chan%s image%s enc%s", ...
%                 num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
%                 num2str(plot_params.image_id), num2str(plot_params.enc_id));
%     mkdir(folder_to_save_in)
%     mkdir(sprintf("%s/single obs/", folder_to_save_in))
%     mkdir(sprintf("%s/mean across trials/", folder_to_save_in))
% 
%     % load PSV
%     PS_file = get_PS_file(params.output_folder, plot_params.patient_id, false);
%     load(PS_file, "label_table")
% 
%     % subset table by plot_params
%     % rows with the image and plot_params
%     if ~without_image
%         rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
%         label_table.patient_ID == plot_params.patient_id & ... % patient
%         label_table.encID_to_imageID(:,plot_params.enc_id)==plot_params.image_id & ... % image in first encoding
%         sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
%     else
%         % rows without the image
%         rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
%         label_table.patient_ID == plot_params.patient_id & ... % patient
%         sum(label_table.encID_to_imageID(:,:)~=plot_params.image_id,2)==3 & ... % image in first encoding
%         sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
%     end
%     subset_table = label_table(rows_to_use,:); % can check
% 
%     % get anat from table
%     plot_params.anat = unique(string(label_table(rows_to_use,:).anatomical_label)); % get anat from table
%     % clear label_table %now would be the time to do it but probably keep for
%     % using trial_id
%     if length(plot_params.anat) > 1
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
%     % flattened_data = possible_PSVs_for_this_obs(plot_params.enc_window_ids, :, :);
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
%                 num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
%                 num2str(plot_params.image_id), num2str(plot_params.enc_id), ...
%                 num2str(subset_table(obs_slice,:).trial_ID));
%         else
%              title_str = sprintf("p%s chan%s WITHOUT image%s enc%s trial%s", ...
%                 num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
%                 num2str(plot_params.image_id), num2str(plot_params.enc_id), ...
%                 num2str(subset_table(obs_slice,:).trial_ID));
%         end
%         fig = plot_PSV(possible_PSVs_for_this_obs(:,:,obs_slice), plot_params, title_str, clims_to_use); % single observation
%         saveas(fig, sprintf("%s/single obs/%s obs%s.png",folder_to_save_in, title_str, num2str(obs_slice)))
%     end
% 
% 
%     % flattened_data = mean(possible_PSVs_for_this_obs(plot_params.enc_window_ids,:,:), 3, 'omitnan');
%     % data_mean = mean(flattened_data(:));
%     % data_std = std(flattened_data(:));
%     % Set the color limits to be 2 standard deviations around the mean
%     % clims_to_use = [data_mean - 2*data_std, data_mean + 2*data_std];
%    % with: [-0.8951    2.1138]; without: [-0.5547    1.3577]
%    clims_to_use = [-1, 2];
% 
%     if ~without_image
%         title_str = sprintf("p%s chan%s image%s enc%s mean across trials", ...
%             num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
%             num2str(plot_params.image_id), num2str(plot_params.enc_id));
%     else
%         title_str = sprintf("p%s chan%s WITHOUT image%s enc%s mean across trials", ...
%             num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
%             num2str(plot_params.image_id), num2str(plot_params.enc_id));
%     end
%     mean(possible_PSVs_for_this_obs, 3, 'omitnan')
%     fig = plot_PSV(mean(possible_PSVs_for_this_obs, 3, 'omitnan'), plot_params, title_str, clims_to_use); % mean of all obs
%     saveas(fig, sprintf("%s/mean across trials/%s.png", folder_to_save_in, title_str))
% 
% end