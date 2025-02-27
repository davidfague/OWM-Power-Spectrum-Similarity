
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