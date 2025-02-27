%% 
% addpath("../subfunctions")
% addpath("../plotting functions")
% addpath("../plotting functions/main")
% addpath("../plotting functions/main nearly")

custom_params = struct();
custom_params.k12wm = false;
custom_params.output_folder_name = 'middle_fixation_baseline';
custom_params.hellbender = true;
custom_params.clip_inf_similarities = true;
params = get_parameters(custom_params);

%% select data
     % patient_id     chan_id      p    anat                      image_id
     % "201915"       42           0    "L Middle Frontal Gyrus"       1 Close p<0.11
     % "201915"       14           0    "L Middle Frontal Gyrus"       2    
     % "201915"       42           0    "L Middle Frontal Gyrus"       3 final cluster good, but raw similarity kinda low
     % "201915"       14           0    "L Middle Frontal Gyrus"       4 
     % "201915"       14           0    "L Middle Frontal Gyrus"       5    
     % "201915"       43           0    "L Middle Frontal Gyrus"       5    
     % "201915"       16           0    "L Middle Frontal Gyrus"       6 Good
     % "201915"       42           0    "L Middle Frontal Gyrus"       6  higher BI similarity   
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
plot_params.patient_id = 8;
plot_params.chan_id = 2;
plot_params.image_id = 8;
plot_params.enc_window_ids = params.enc1_win_IDs; % change with enc_id
plot_params.enc_id = 1;
params.session_id = 1;

% for PSVs
plot_params.frequencies_to_use = 1:40; % empty for whatever the data is % for plotting PSVs only

% for WI vs BI
plot_params.type = 'EMS'; % for plotting WI_vs_BI only
plot_params.mean_out_time_dimensions = false;
plot_params.average_diff = true;
plot_params.same_n = true; % only affects plot_similarity_means_heatmaps

plot_params.only_all3_correct = true; % filters test trials (with item)

% es_freq_bands = {[1:8], [8:20], [20:40], [1:40]};
es_freq_bands = {[30:70]};
% es_freq_bands = {[1:40]};
% look at power spectra to determine bands
% 4:8
% 8:12-14
% 15:30
% 30:70
% 70:140

% mean spectra over time (whole trial(after fix)) and create some figures of different areas
% mean,std across patients, channels for each area % mtg, mfg, hippocampus,
% amygdala

% correct plot
%%
plot_params.recompute_table = true;
% params.patient_IDs = [010];

if plot_params.recompute_table
    original_mean_out_time = plot_params.mean_out_time_dimensions;
    plot_params.mean_out_time_dimensions = true;
    for freqs_idx = 1:length(es_freq_bands)
        all_p_table = table();
        params.ES_freq_band = es_freq_bands{freqs_idx};
        params.freq_min = min(params.ES_freq_band);
        params.freq_max = max(params.ES_freq_band);

        for pat_idx = 1:length(params.patient_IDs)
            plot_params.patient_id = params.patient_IDs(pat_idx);
            all_p_table = calc_WI_vs_BI_table(params, plot_params, all_p_table);
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
plot_params.mean_out_time_dimensions = original_mean_out_time;
clearvars patient_ids_to_use table_all_patients original_mean_out_time

%% compute
% plot_params.patient_id = 201908;
% plot_params.image_id = 7;
% plot_params.chan_id = 17;
% for row_idx = 1:length(data.patient_id)
    % plot_params.chan_id = data(row_idx,:).chan_id;
    % plot_params.image_id = data(row_idx,:).image_id;
    % plot_params.patient_id = data(row_idx,:).patient_id;
    % fprintf("row%s p%s chan%s image%s\n",num2str(row_idx), num2str(plot_params.patient_id), num2str(plot_params.chan_id), num2str(plot_params.image_id))
    % % plot PSVs
    % WI
    plot_PSVs(params, plot_params, false) % plot trials with the image
    close all
    
    % BI
    plot_PSVs(params, plot_params, true) % plot trials without the image
    close all
    
    % % plot WI vs BI EMS, t-test, clusters
    % plot WI vs BI EMSs
    % for ave_diff = [true, false]
    
        % plot_params.average_diff = ave_diff;
    [real_t_values, real_p_values, surrogate_t_values, surrogate_p_values] = plot_WI_vs_BI(params, plot_params);
    close all
% end

%% plot within_trial
% plot_params.patient_id = 201908;
% plot_params.image_id = 7;
% plot_params.chan_id = 17;
% plot_params.without_image = false;
% plot_within_trial(params, plot_params)
% close all

% for row_idx = 8:length(data.patient_id)
    % plot_params.chan_id = data(row_idx,:).chan_id;
    % plot_params.image_id = data(row_idx,:).image_id;
    % plot_params.patient_id = data(row_idx,:).patient_id;

    % with image
    for without_image = [true, false]
        plot_params.without_image = without_image;
        plot_within_trial(params, plot_params)
        close all
    end

% end

%% functions

function all_p_table = calc_WI_vs_BI_table(params, plot_params, all_p_table)
    patient_id = [plot_params.patient_id];

        % load patient's preprocessed data
    % patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
    patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, patient_id);

    for session_idx = 1:length(patient_preprocessed_data_paths)
        params.session_id = session_idx;
        patient_preprocessed_data_path = patient_preprocessed_data_paths{session_idx};
        disp(patient_preprocessed_data_path)

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
    
        for image_id = 1:9
            plot_params.image_id = image_id;
            channel_ids_to_use = extractIntegersFromFilenames(patient_id, params.comp_options, plot_params.enc_id, image_id, params); % image_id=1 is used as a default here. Channels should change with respect to iamge_id.
    
            for chan_idx = 1:length(channel_ids_to_use)
                chan_id = channel_ids_to_use(chan_idx);
                plot_param.chan_id = chan_id;
                fprintf("p%s session%d chan%s image%s %d-%dHz\n", num2str(patient_id), session_idx, num2str(chan_id), num2str(image_id), params.freq_min, params.freq_max);
                anat = string(anat_labels.labelsanatbkedit.anatmacro1(chan_id));

                % get WI, BI matrices; size: (nEtimes, nMtimes, nTrialCombinations)
                [WI, BI] = get_and_pre_process_WI_BI(patient_id, params.comp_options, plot_params.enc_id, ...
                    image_id, chan_id, params, plot_params.type, plot_params.only_all3_correct);
    
                if numel(BI) == 0 || numel(WI) == 0
                    fprintf("issue with this BI. Likely some some nans everywhere... skipping\n")
                    final_p = nan;
                    result_table = table(patient_id, ...
                                     chan_id, final_p, anat, ...
                                     image_id, params.session_id,...
                                     'VariableNames', {'patient_id', 'chan_id', 'p', 'anat', 'image_id', 'session_id'});
                
                    all_p_table = [all_p_table; result_table];
                    continue
                end
        
                if params.clip_inf_similarities
                   [WI, BI] = clip_infs_of_z_similarities(WI,BI);
                end
        
                % temporal generalization
                if plot_params.mean_out_time_dimensions
                    WI=mean(WI, [1 2]);
                    BI=mean(BI, [1 2]);
                end
        
                % calculate real differences
                fprintf("    calculating real differences\n")
                if plot_params.average_diff
                    diff = calc_diff_avg(WI, BI, 1000);
                else
                    diff = calc_diff_all_pairs(WI, BI, 10000);
                end
        
                fprintf("    calculating real t-tests\n")
                [real_t_values, real_p_values] = calc_ttest(diff);
            
                % calculate 1000 surrogate differences and t-tests
                fprintf("    calculating surrogate differences and t-tests\n")
                [surrogate_t_values, surrogate_p_values, surrogate_diff] = calc_surrogate(WI, BI, 1000, plot_params.average_diff, plot_params.mean_out_time_dimensions);
                % save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in),"real_p_values", "real_t_values", "surrogate_t_values", "surrogate_p_values", "surrogate_diff")
        
                if any(isnan(real_t_values))
                    error("nans detected in real_t_values.")
                elseif any(isnan(real_p_values))
                    error("nans detecte in real_p_values")
                elseif any(isnan(surrogate_p_values))
                    error("nans detected in surrogate_p_values")
                elseif any(isnan(surrogate_t_values))
                    error("nans detected in surrogate_t_values")
                end

                final_p = mean(real_t_values > surrogate_t_values);
        
                    % Create the table
                result_table = table(patient_id, ...
                                     chan_id, final_p, anat, ...
                                     image_id, ...
                                     session_idx, ...
                                     'VariableNames', {'patient_id', 'chan_id', 'p', 'anat', 'image_id', 'session_id'});
                
                all_p_table = [all_p_table; result_table];
            end
        end
    end
end

function [real_t_values, real_p_values, surrogate_t_values, surrogate_p_values, final_p] = plot_WI_vs_BI(params, plot_params)
    tic
    patient_id = [plot_params.patient_id];
    channel_ids_to_use = [plot_params.chan_id];
    patient_preprocessed_data_paths = get_patient_preprocessed_data_path(params, patient_id);

    for session_idx = 1:length(patient_preprocessed_data_paths)
        params.session_id = session_idx;
        patient_preprocessed_data_path = patient_preprocessed_data_paths{session_idx};
        
        plot_params.WI_BI_folder_to_save_in = sprintf("results/WI vs BI/p%s chan%s image%s enc%s sess%d", ...
                    num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
                    num2str(plot_params.image_id), num2str(plot_params.enc_id), session_idx);
    
        folder_to_save_in = plot_params.WI_BI_folder_to_save_in;
    
        warning('off', 'MATLAB:MKDIR:DirectoryExists');
        mkdir(folder_to_save_in);
        warning('on', 'MATLAB:MKDIR:DirectoryExists');
    
    
        % patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
        % image_labels = load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C'); % if wanted
        if params.k12wm
            prefix = strsplit(patient_preprocessed_data_path,'\');
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
        
        for chan_idx = 1:length(channel_ids_to_use)
            chan_id = channel_ids_to_use(chan_idx);
            plot_params.anat = string(anat_labels.labelsanatbkedit.anatmacro1(chan_id));
            fprintf("   %s\n", plot_params.anat)
    
            % get WI, BI matrices; size: (nEtimes, nMtimes, nTrialCombinations)
            [WI, BI] = get_and_pre_process_WI_BI(plot_params.patient_id, params.comp_options, plot_params.enc_id, ...
                plot_params.image_id, plot_params.chan_id, params, plot_params.type, plot_params.only_all3_correct);
    
            % load some previous p vlaues if needed @DEPRECATING
            % selected_p_data = load(sprintf("results/selected_p_data_p%s_image%s_chan%s.mat", ...
            %     num2str(plot_params.patient_id),num2str(plot_params.image_id),num2str(plot_params.chan_id)));
    
            % plot EMS timecourse
            fig = plot_similarity_means_timecourse(WI,BI, plot_params);
            saveas(fig, sprintf("%s/WI_BI_time_courses_same_n.fig", folder_to_save_in))
    
            % plot EMS heatmap
            for same_n = [true, false]
                plot_params.same_n = same_n;
                fig = plot_similarity_means_heatmap(WI,BI, plot_params);
                if plot_params.same_n
                    saveas(fig, sprintf("%s/WI_BI_heatmaps_same_n.fig", folder_to_save_in))
                else
                    saveas(fig, sprintf("%s/WI_BI_heatmaps.fig", folder_to_save_in))
                end
            end
    
            close all % close first 3 figures
    
            % temporal generalization
            if plot_params.mean_out_time_dimensions
                WI=mean(WI, [1 2]);
                BI=mean(BI, [1 2]);
            end
    
            % calculate real differences
            fprintf("    calculating real differences\n")
            if plot_params.average_diff
                diff = calc_diff_avg(WI, BI, 1000);
            else
                diff = calc_diff_all_pairs(WI, BI, 10000);
            end
    
            if ~exist(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in), 'file')
                % perform t-test on real differences
                fprintf("    calculating real t-tests\n")
                [real_t_values, real_p_values] = calc_ttest(diff);
            
                % calculate 1000 surrogate differences and t-tests
                fprintf("    calculating surrogate differences and t-tests\n")
                [surrogate_t_values, surrogate_p_values, surrogate_diff] = calc_surrogate(WI, BI, 1000, plot_params.average_diff, plot_params.mean_out_time_dimensions);
                save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in),"real_p_values", "real_t_values", "surrogate_t_values", "surrogate_p_values", "surrogate_diff")
            else
                load(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in))
                % if ~exist("surrogate_diff","var")
                %     [surrogate_t_values, surrogate_p_values, surrogate_diff] = calc_surrogate(WI, BI, 1000, plot_params.average_diff, plot_params.mean_out_time_dimensions);
                %     save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in),"real_p_values", "real_t_values", "surrogate_t_values", "surrogate_p_values", "surrogate_diff")
                % end
            end
               
            % % plot t_values
            % fig = plot_t_values(real_value, null_values);
    
            % perform clustering and get final p or get final p from temporal
            % generalization
            if ~exist("final_p","var") || ~exist("all_real_clusters","var")
                if ~plot_params.mean_out_time_dimensions
                    % clustering p
                    [final_p, all_real_clusters, significant_clusters, sum_threshold, alpha] = cluster_analysis_and_visualization_abs_main(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values, plot_params);
                    save(sprintf("%s/significant_cluster_data.mat", plot_params.WI_BI_folder_to_save_in), "-append", "significant_clusters", "sum_threshold", "alpha", "final_p", "all_real_clusters")
                else
                    % temproal generalization
                    % final_p = mean(real_p_values > surrogate_p_values);
                    final_p = mean(real_t_values > surrogate_t_values);
                end
            end
    
            % plot significant clusters
            fig = plot_WI_BI_significant_clusters(BI,WI, significant_clusters, all_real_clusters, plot_params);
            saveas(fig, sprintf("%s/main_cluster_visual.fig", folder_to_save_in))
            saveas(fig, sprintf("%s/main_cluster_visual.png", folder_to_save_in))
    
            % plot differences distributions
            % if ~plot_params.mean_out_time_dimensions
            %     fig = plot_diff(diff, surrogate_diff(:), final_p, plot_params);
            %     if plot_params.average_diff
            %         saveas(fig, sprintf("%s/WI_minus_BI_distributions_trialavgs.fig", folder_to_save_in))
            %     else
            %         saveas(fig, sprintf("%s/WI_minus_BI_distributions_trial_individuals.fig", folder_to_save_in))
            %     end
            % end
    
        end
        fprintf("finished computing in %.2f seconds\n",toc);
        close all
    end
end

function plot_within_trial(params, plot_params)
    PS_file = get_PS_file(params, plot_params.patient_id, false);
    load(PS_file, 'label_table');
    ES_file = get_ES_file(PS_file, false, false); % dims: rows=10440,          encT=31         trialT=630           encsN=3
    plot_params.WTS_folder_to_save_in = sprintf("results/WI vs BI/p%s chan%s image%s enc%s/WTS/", ...
            num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
            num2str(plot_params.image_id), num2str(plot_params.enc_id));

    fprintf("plotting WTS for %s\n", plot_params.WTS_folder_to_save_in)

    folder_to_save_in = plot_params.WTS_folder_to_save_in;
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir(folder_to_save_in);
    warning('on', 'MATLAB:MKDIR:DirectoryExists');

    if ~plot_params.without_image
        rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
        label_table.patient_ID == str2double(plot_params.patient_id) & ... % patient
        label_table.encID_to_imageID(:,plot_params.enc_id)==plot_params.image_id & ... % image in first encoding
        sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
    else
        % rows without the image
        rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
        label_table.patient_ID == str2double(plot_params.patient_id) & ... % patient
        sum(label_table.encID_to_imageID(:,:)~=plot_params.image_id,2)==3 & ... % image in first encoding
        sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
    end

    subset_table = label_table(rows_to_use,:); % can check
    plot_params.anat = string(unique(string(subset_table.anatomical_label)));

    load(ES_file, 'all3_ES_matrix')
    ES_matrix = all3_ES_matrix(rows_to_use,:,:,plot_params.enc_id);

    rows_without_nans = ~all(isnan(ES_matrix),[2 3]);
    ES_matrix = ES_matrix(rows_without_nans,:,:); % remove NANs
    subset_table = subset_table(rows_without_nans,:);
    clear all3_ES_matrix

    max_obs_to_plot = 3;
    n_obs_to_plot = min(max_obs_to_plot, size(subset_table, 1));

    % plot individual observations
    if plot_params.without_image
        with_image_str = 'WITHOUT';
    else
        with_image_str = 'WITH';
    end

    enforced_clim = [-0.3 0.6];

    for row_idx=1:n_obs_to_plot
        plot_params.trial_id = subset_table.trial_ID(row_idx);
        title_str = sprintf("p%s chan%s %s image%s TRIAL%s", num2str(plot_params.patient_id), ...
        num2str(plot_params.chan_id), with_image_str, num2str(plot_params.image_id), ...
        num2str(plot_params.trial_id));
        fig = plot_within_trial_similarity(ES_matrix(row_idx,:,:), plot_params, title_str, enforced_clim);
        savefig(fig, sprintf("%s/TRIAL %s %s.fig", folder_to_save_in, num2str(plot_params.trial_id), with_image_str))
    end

    % plot mean obs
    plot_params.trial_id = nan;
    title_str = sprintf("p%s chan%s %s image%s MEAN ACROSS TRIALS n=%s", num2str(plot_params.patient_id), ...
        num2str(plot_params.chan_id), with_image_str, num2str(plot_params.image_id), num2str(size(ES_matrix,1)));
    fig = plot_within_trial_similarity(mean(ES_matrix,1), plot_params, title_str, enforced_clim);
    savefig(fig, sprintf("%s/MEAN ACROSS TRIALS %s.fig", folder_to_save_in, with_image_str))
    saveas(fig, sprintf("%s/MEAN ACROSS TRIALS %s.png", folder_to_save_in, with_image_str))

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