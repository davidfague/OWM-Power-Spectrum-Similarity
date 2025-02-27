function plot_PSVs(params, plot_params, without_image)
% plots the whole-trial PSVS for the individual trial or mean across trials for a specific chan, image,
% enc slice.
    folder_to_save_in = sprintf("results/PSV plots/same_clims/p%s chan%s image%s enc%s", ...
                num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
                num2str(plot_params.image_id), num2str(plot_params.enc_id));
    mkdir(folder_to_save_in)
    mkdir(sprintf("%s/single obs/", folder_to_save_in))
    mkdir(sprintf("%s/mean across trials/", folder_to_save_in))
    
    % load PSV
    PS_file = get_PS_file(params, plot_params.patient_id, false);
    load(PS_file, "label_table")
    
    % subset table by plot_params
    if ~without_image 
        % rows with the image and plot_params
        rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
        label_table.patient_ID == plot_params.patient_id & ... % patient
        label_table.encID_to_imageID(:,plot_params.enc_id)==plot_params.image_id & ... % image in first encoding
        sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
    else
        % rows without the image
        rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
        label_table.patient_ID == plot_params.patient_id & ... % patient
        sum(label_table.encID_to_imageID(:,:)~=plot_params.image_id,2)==3 & ... % image in first encoding
        sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
    end
    % rows_to_use = ones(size(rows_to_use));

    subset_table = label_table(rows_to_use,:); % can check
    
    % get anat from table
    plot_params.anat = unique(string(label_table(rows_to_use,:).anatomical_label)); % get anat from table
    % clear label_table %now would be the time to do it but probably keep for
    % using trial_id
    if length(plot_params.anat) > 1
        error("should only detect one anat")
    end
    
    % get correspoding rows (trials) of PSVs from PS_file
    % PS_file = matfile(PS_file);
    load(PS_file, 'all_windowed_mean_PS_vectors')
    possible_PSVs_for_this_obs = all_windowed_mean_PS_vectors(:,:,rows_to_use);
    rows_without_nans =~any(isnan(possible_PSVs_for_this_obs),[1 2]);
    possible_PSVs_for_this_obs = possible_PSVs_for_this_obs(:,:,rows_without_nans); % remove trials with nans
    subset_table = subset_table(rows_without_nans,:);
    clear all_windowed_mean_PS_vectors
    
    flattened_data = possible_PSVs_for_this_obs(plot_params.enc_window_ids, :, :);
    % Calculate the mean and standard deviation of the flattened data
    data_mean = mean(flattened_data(:));
    data_std = std(flattened_data(:));
    % Set the color limits to be 2 standard deviations around the mean
    clims_to_use = [data_mean - 2*data_std, data_mean + 2*data_std];
    % clims_to_use = [-1, 2];% with:[-3.3583    4.5770]; without: [-3.5836    4.3866]
    
    max_iter = 9;
    max_iter_to_use = min(max_iter, size(possible_PSVs_for_this_obs,3));
    for obs_slice = 1:max_iter_to_use % can also select 1
        if ~without_image
            title_str = sprintf("p%s chan%s image%s enc%s trial%s", ...
                num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
                num2str(plot_params.image_id), num2str(plot_params.enc_id), ...
                num2str(subset_table(obs_slice,:).trial_ID));
        else
             title_str = sprintf("p%s chan%s WITHOUT image%s enc%s trial%s", ...
                num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
                num2str(plot_params.image_id), num2str(plot_params.enc_id), ...
                num2str(subset_table(obs_slice,:).trial_ID));
        end
        fig = plot_PSV(possible_PSVs_for_this_obs(:,:,obs_slice), plot_params, title_str, clims_to_use); % single observation
        saveas(fig, sprintf("%s/single obs/%s obs%s.png",folder_to_save_in, title_str, num2str(obs_slice)))
    end
    
    
    % flattened_data = mean(possible_PSVs_for_this_obs(plot_params.enc_window_ids,:,:), 3, 'omitnan');
    % data_mean = mean(flattened_data(:));
    % data_std = std(flattened_data(:));
    % Set the color limits to be 2 standard deviations around the mean
    % clims_to_use = [data_mean - 2*data_std, data_mean + 2*data_std];
   % with: [-0.8951    2.1138]; without: [-0.5547    1.3577]
   clims_to_use = [-1, 2];
    
    if ~without_image
        title_str = sprintf("p%s chan%s image%s enc%s mean across trials", ...
            num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
            num2str(plot_params.image_id), num2str(plot_params.enc_id));
    else
        title_str = sprintf("p%s chan%s WITHOUT image%s enc%s mean across trials", ...
            num2str(plot_params.patient_id), num2str(plot_params.chan_id), ...
            num2str(plot_params.image_id), num2str(plot_params.enc_id));
    end
    mean(possible_PSVs_for_this_obs, 3, 'omitnan')
    fig = plot_PSV(mean(possible_PSVs_for_this_obs, 3, 'omitnan'), plot_params, title_str, clims_to_use); % mean of all obs
    saveas(fig, sprintf("%s/mean across trials/%s.png", folder_to_save_in, title_str))

end