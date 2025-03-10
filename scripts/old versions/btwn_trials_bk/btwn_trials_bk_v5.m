%% between trials comparisons - within-item and between-items
% note: 'test trials' refer to trials that had encoding taken; 'control trials' refers to trials that had maintenance taken

% v5 changing plotting: look at individual trial pairs for the same 'test'
% trial, create table to summarize plots across all anat and all test trials;
% look at within-image, between-trial similarity as function of trial_id ('real-time')
% v4 update file names, change some plotting
% v3 enabling params; remove nan filtering based on table.
% v2 changing it so that BT-WI and BT-BI are subplots 2 and 3
% v2 additional: adding sig. cluster analysis between BT-WI and BT-BI
% v1 WT-WI is subplot1, BT-WI or BT-BI is subplot 2.

params = get_parameters();

% Set patient, image, and encoding IDs
patient_id = 201901;
image_id = 1;
enc_id = 1;
channel_ids_to_use = [10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66];

% Comparison options and titles
comp_options = {'corr BT BI', 'corr BT WI'};
comp_titles = {'Between Trials - Between Images', 'Between Trials - Within Images'};

old=true; % whether to use true:(100 ms window, -.25--.75 baseline within trial) or false:(200 ms, -0.5-0.0 baseline across trials)
if old
    output_folder_name = 'allpatients gammamod allregions allitem allenc';
    params.output_destination = "C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBTvsBI newer_baselining_windowing"; % Set your desired output folder
else
    output_folder_name = 'allpatients gammamod allregions allitem allenc baseline across trials';
    params.output_destination = 'C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBTvsBI original_baselining_windowing'; % Set your desired output folder
end

params.output_destination = fullfile(params.output_destination, sprintf("p%s_enc%s_image%s", num2str(patient_id), num2str(enc_id), num2str(image_id)));
clear old

% Load required data
load(fullfile(strcat(params.output_folder,"\\",num2str(patient_id),"\\encoding_similarity.mat"))) % within trial similarity
load(fullfile(strcat(params.output_folder,"\\",num2str(patient_id),"\\PSVs.mat")), 'label_table') % label table
patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C');
clear patient_preprocessed_data_path

% rows_without_nan_by_enc = ~squeeze(any(any(isnan(all3_ES_matrix), 2), 3));

%% identify Unique anatomical labels
[unique_channels_labels, ia] = unique(label_table.channel_ID, 'rows');
% rows_without_nan = get_rows_without_nan(label_table);
% label_table = label_table(rows_without_nan, :);

%% identify Channels to use
% channels_to_use = contains(string(label_table{ia, 3}), 'amyg', 'IgnoreCase', true) | ...
%                   contains(string(label_table{ia, 3}), 'hipp', 'IgnoreCase', true);
% channels_to_use = label_table{ia(channels_to_use), 2};

% amyg, hipp channels
% for chan_id = [10 11 12  22  61] % 201901
% for chan_id = [1 2  10 11 12 13  46 47 48 49 50  57 58 66] % 201910

% mfg, mtg channels
% for chan_id = [7 8 9   22 23 24 25 26   48 49 50   63 64 65 66] % 201901
% for chan_id = [7 14  25 26 27  41 42  52 53 54  59 60] % 201910
% for chan_id = [90 18 26 28 36 39] % was for 201907 but weird data

%% plot single-trial WT vs BTWI vs BTBI heatmaps
[all_channels_data] = load_all_channels_data(params, patient_id, comp_options, enc_id, image_id);

for chan_id = channel_ids_to_use

    % Load between_trials info for the channel (should be same whether
    % loaded from comparison{1} or comparison{2}

    % Load BT-BI data for this channel
    [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);

    %size(BT_ES) = 31   380    90    12 (enc, maint, ntrials_control, ntrials_test)


    Mtime = [250:630]; % Delay time

    figure;

    sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
        patient_id, chan_id, enc_id, C{1, image_id}), ...
        unique(string(BT_WI_data.chan_test_table.anatomical_label))]);

    % get test trial
    test_trial_to_use = 1; % first, second, etc.; change if nan % test refers to target_item_correct_trials
    original_row_for_test_trial_1 = find((label_table.trial_ID == BT_BI_data.chan_test_table(test_trial_to_use,:).trial_ID) & (label_table.channel_ID == BT_BI_data.chan_test_table(test_trial_to_use,:).channel_ID));
    within_trial_single_trial = squeeze(all3_ES_matrix(original_row_for_test_trial_1, :, Mtime, enc_id)); % result = all3_ES_matrix(rows_without_nan, :, Mtime, enc_id);
    if any(isnan(within_trial_single_trial(:)))
        error("nan in this within trial data. change test_trial_to_use")
    end

    % define control trials to use
    control_trials_to_use = 2:size(BT_WI_data.chan_test_table, 1);

    % Subplot 1: Within-trial similarity for test trial
    subplot(1+length(control_trials_to_use), 2, 1)
    imagesc(within_trial_single_trial, [0 0.5])
    title('Within-Trial Similarity - single trial')
    add_single_trial_WT_vs_BTWI_vs_BTBI_labels()

    % subplots 2-3 and on, BT data
    iter=1; % make subplotting easier
    for control_trial_to_use = control_trials_to_use
        if test_trial_to_use == control_trial_to_use
            error("test_trial_to_use and second_test_trial_to_use_for_pair are the same. BT-WI will become WT-WI (misleading)")
        end

        % Subplot 2: Between-trial Within-Image single-trial pair
        subplot(1+length(control_trials_to_use), 2, 2+iter)
        imagesc(squeeze(BT_WI_data.BT_ES(:,:,control_trial_to_use,test_trial_to_use)), [0 0.5])
        title(strcat(comp_titles{2}, ' same trial, 1-trial pairs n=', num2str(size(BT_WI_data.BT_ES(:,:,control_trial_to_use,test_trial_to_use),3))))
        add_single_trial_WT_vs_BTWI_vs_BTBI_labels()
    
        % Subplot 3: Between-trial Between-Image single-trial pair 
        subplot(1+length(control_trials_to_use), 2, 3+iter)
        imagesc(squeeze(BT_BI_data.BT_ES(:,:,control_trial_to_use, test_trial_to_use)), [0 0.5])
        title(strcat(comp_titles{1}, ' same trial, 1-trial pairs n=', num2str(size(BT_BI_data.BT_ES(:,:,control_trial_to_use,test_trial_to_use),3))))
        add_single_trial_WT_vs_BTWI_vs_BTBI_labels()
        % BT_BI_data.chan_ctrl_table
        iter=iter+2;
    end
    save_file = fullfile(params.output_destination, sprintf('single_trial_WT_vs_BTWI_vs_BTBI_heatmaps/Patient%d_Channel%d_Enc%d_Image%d.fig', ...
    patient_id, chan_id, enc_id, image_id));

    save_dir = fileparts(save_file);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    saveas(gcf, save_file);

end
function add_single_trial_WT_vs_BTWI_vs_BTBI_labels()
    % ylabel('Encoding Window ID')
    % xlabel('Maintenance Window ID')
    colorbar;

end

%% plot single-trial WT vs BTWI vs BTBI lines (meaning across encoding)
same_n = true;
[all_channels_data] = load_all_channels_data(params, patient_id, comp_options, enc_id, image_id);

% summary_table = table();

for chan_id = channel_ids_to_use

    % Load between_trials info for the channel (should be same whether
    % loaded from comparison{1} or comparison{2}

    % Load BT-BI data for this channel
    [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);
    %size(BT_ES) = 31   380    90    12 (enc, maint, ntrials_control, ntrials_test)

    flattenedMatrix = reshape(BT_BI_data.BT_ES, [], size(BT_BI_data.BT_ES, 3), size(BT_BI_data.BT_ES, 4));
    flattenedMatrix = squeeze(any(isnan(flattenedMatrix),1));
    test_trials_without_nan = ~all(flattenedMatrix, 1);
    ctrl_trials_without_nan = ~all(flattenedMatrix, 2);
    % figure;
    % imagesc(flattenedMatrix)

    % rows_without_nan = ~squeeze(any(isnan(BT_BI_data.BT_ES), [1,2,4]));
    BT_BI_data.BT_ES = BT_BI_data.BT_ES(:,:,ctrl_trials_without_nan,test_trials_without_nan);
    % rows_without_nan = ~squeeze(any(isnan(BT_WI_data.BT_ES), [1,2,4]));
    BT_WI_data.BT_ES = BT_WI_data.BT_ES(:,:,test_trials_without_nan,test_trials_without_nan);



    Mtime = [250:630]; % Delay time

    figure;
    hold on

    sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
        patient_id, chan_id, enc_id, C{1, image_id}), ...
        unique(string(BT_WI_data.chan_test_table.anatomical_label))]);

    % get test trial
    test_trial_to_use = 1; % first, second, etc.; change if nan % test refers to target_item_correct_trials
    % test_trials_to_use = 1:size(BT_BI_data.BT_ES, 4); % update to loop through
    original_row_for_test_trial_1 = find((label_table.trial_ID == BT_BI_data.chan_test_table(test_trial_to_use,:).trial_ID) & (label_table.channel_ID == BT_BI_data.chan_test_table(test_trial_to_use,:).channel_ID));
    within_trial_single_trial = squeeze(all3_ES_matrix(original_row_for_test_trial_1, :, Mtime, enc_id)); % result = all3_ES_matrix(rows_without_nan, :, Mtime, enc_id);
    if any(isnan(within_trial_single_trial(:)))
        error("nan in this within trial data. change test_trial_to_use")
    end

    % define control trials to use
    % control_trials_to_use = 2:size(BT_WI_data.chan_test_table, 1);
    % control_trials_to_use_BT_BI = 1:size(BT_BI_data.chan_ctrl_table, 1);
    control_trials_to_use = 2:size(BT_WI_data.BT_ES, 3);
    control_trials_to_use_BT_BI = 1:size(BT_BI_data.BT_ES, 3);

    % Subplot 1: Within-trial similarity for test trial
    subplot(1, 3, 1)
    plot(squeeze(mean(within_trial_single_trial,1)))
    title('Within-Trial Similarity - single trial')
    add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
    hold on

    % subplots 2-3 and on, BT data
    for control_trial_to_use = control_trials_to_use
        if test_trial_to_use == control_trial_to_use
            error("test_trial_to_use and second_test_trial_to_use_for_pair are the same. BT-WI will become WT-WI (misleading)")
        end

        % Subplot 2: Between-trial Within-Image single-trial pair
        subplot(1, 3, 2)
        plot(squeeze(mean(BT_WI_data.BT_ES(:,:,control_trial_to_use,test_trial_to_use), 1)))
        % title(strcat(comp_titles{2}, ' same trial, 1-trial pairs n=', num2str(size(BT_WI_data.BT_ES(:,:,control_trial_to_use,test_trial_to_use),3))))
        add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
        hold on
    
        if same_n
            % Subplot 3: Between-trial Between-Image single-trial pair 
            subplot(1, 3, 3)
            plot(squeeze(mean(BT_BI_data.BT_ES(:,:,control_trial_to_use, test_trial_to_use), 1)))
            % title(strcat(comp_titles{1}, ' same trial, 1-trial pairs n=', num2str(size(BT_BI_data.BT_ES(:,:,control_trial_to_use,test_trial_to_use),3))))
            add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
            % BT_BI_data.chan_ctrl_table
            hold on
        end
    end
    if ~same_n
            subplot(1, 3, 3)
            plot(squeeze(mean(BT_BI_data.BT_ES(:,:,control_trials_to_use_BT_BI, test_trial_to_use), 1)))
            % title(strcat(comp_titles{1}, ' same trial, 1-trial pairs n=', num2str(size(BT_BI_data.BT_ES(:,:,control_trialw_to_use_BT_BI,test_trial_to_use),3))))
            add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
            % BT_BI_data.chan_ctrl_table
            hold on
    end

    % plot horizontal line for the average
        % Subplot 1: Within-trial similarity for test trial
        subplot(1, 3, 1)
        WT_mean = squeeze(mean(within_trial_single_trial,[1, 2]));
        yline(WT_mean)
        title(strcat('Within-Trial Similarity - single trial, mean: ', num2str(WT_mean), ' n:1'))
        add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
        hold on
        % Subplot 2: Between-trial Within-Image single-trial pair
        subplot(1, 3, 2)
        BT_WI_mean = squeeze(mean(BT_WI_data.BT_ES(:,:,control_trials_to_use,test_trial_to_use), [1,2,3]));
        yline(BT_WI_mean)
        n_BT_WI = length(control_trials_to_use);
        title(strcat(comp_titles{2}, ' mean:  ', num2str(BT_WI_mean), ' n: ', num2str(n_BT_WI)))
        add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
        hold on
        % Subplot 3: Between-trial Between-Image single-trial pair 
        subplot(1, 3, 3)
        if same_n
            BT_BI_mean = squeeze(mean(BT_BI_data.BT_ES(:,:,control_trials_to_use, test_trial_to_use), [1,2,3]));
            n_BT_BI = length(control_trials_to_use);
        else
            BT_BI_mean = squeeze(mean(BT_BI_data.BT_ES(:,:,:, test_trial_to_use), [1,2,3]));
            n_BT_BI = length(control_trials_to_use_BT_BI);
        end
        title(strcat(comp_titles{1}, ' mean:  ', num2str(BT_BI_mean), ' n: ', num2str(n_BT_BI)))
        yline(BT_BI_mean)
        add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
        hold on

    % add to table
    % anat = unique(string(BT_WI_data.chan_test_table.anatomical_label));
    % new_row = table(anat, WT_mean, BT_WI_mean, BT_BI_mean, n_BT_WI, n_BT_BI, ...
    %             'VariableNames', {'anat', 'mean_WT', 'mean_BT_WI', 'mean_BT_BI', 'n_BT_WI', 'n_BT_BI'});
    % summary_table = [summary_table; new_row]; % Append the new row

    if ~exist(params.output_destination, 'dir')
        mkdir(params.output_destination)
    end
    if same_n
        save_file = fullfile(params.output_destination, sprintf('single_trial_WT_vs_BTWI_vs_BTBI_lines_same_n/Patient%d_Channel%d_Enc%d_Image%d.fig', ...
        patient_id, chan_id, enc_id, image_id));
    else
        save_file = fullfile(params.output_destination, sprintf('single_trial_WT_vs_BTWI_vs_BTBI_lines/Patient%d_Channel%d_Enc%d_Image%d.fig', ...
        patient_id, chan_id, enc_id, image_id));
    end

    save_dir = fileparts(save_file);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    saveas(gcf, save_file);
end

% summary_table.WIvsBI = summary_table.mean_BT_WI - summary_table.mean_BT_BI;
% save(fullfile(params.output_destination, "WIvsBI_summary_table.mat"), "summary_table")

function add_single_trial_WT_vs_BTWI_vs_BTBI_labels_line()
    % ylabel('Encoding Window ID')
    % xlabel('Maintenance Window ID')
    ylim([0, 1])

end

function [all_channels_data] = load_all_channels_data(params, patient_id, comp_options, enc_id, image_id)
    all_channels_data = load(sprintf('%s\\%s\\%s\\enc%d_image%d\\all_channels_data.mat', ...
        params.output_folder, num2str(patient_id), comp_options{2}, enc_id, image_id));
end

function [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params)
    % Load BT-BI data for this channel
    BT_BI_data = load(sprintf('%s\\%s\\%s\\enc%s_image%s\\BT_%s.mat', ...
        params.output_folder, num2str(patient_id), comp_options{1}, num2str(enc_id), num2str(image_id), num2str(chan_id)));
    clear BT_ES
    % BT_BI_data_std = control_EMS_matrices_std; % Rename for clarity

    % Load BT-WI data for this channel
    BT_WI_data = load(sprintf('%s\\%s\\%s\\enc%s_image%s\\BT_%s.mat', ...
        params.output_folder, num2str(patient_id), comp_options{2}, num2str(enc_id), num2str(image_id), num2str(chan_id)));
end

%% only table no plotting

summary_table = table();

for chan_id = channel_ids_to_use
    % Load BT-BI and BT-WI data for the channel
    [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);

    % Preprocess data to remove NaNs
    flattenedMatrix = reshape(BT_BI_data.BT_ES, [], size(BT_BI_data.BT_ES, 3), size(BT_BI_data.BT_ES, 4));
    flattenedMatrix = squeeze(any(isnan(flattenedMatrix), 1));
    test_trials_without_nan = ~all(flattenedMatrix, 1);
    ctrl_trials_without_nan = ~all(flattenedMatrix, 2);
    BT_BI_data.BT_ES = BT_BI_data.BT_ES(:, :, ctrl_trials_without_nan, test_trials_without_nan);
    BT_WI_data.BT_ES = BT_WI_data.BT_ES(:, :, test_trials_without_nan, test_trials_without_nan);

    test_trials_to_use = 1:size(BT_BI_data.BT_ES, 4); % Loop through all test trials

    % Iterate over test trials
    for test_trial_to_use = test_trials_to_use
        % Skip if the test trial has NaNs
        % if any(isnan(squeeze(BT_WI_data.BT_ES(:, :, :, test_trial_to_use)), 'all'))
        %     continue;
        % end

        % Compute means for the current test trial
        control_trials_to_use = 2:size(BT_WI_data.BT_ES, 3);
        control_trials_to_use_BT_BI = 1:size(BT_BI_data.BT_ES, 3);

        WT_mean = squeeze(mean(BT_WI_data.BT_ES(:, :, test_trial_to_use, test_trial_to_use), [1, 2]));
        BT_WI_mean = squeeze(mean(BT_WI_data.BT_ES(:, :, control_trials_to_use, test_trial_to_use), [1, 2, 3]));

        if same_n
            BT_BI_mean = squeeze(mean(BT_BI_data.BT_ES(:, :, control_trials_to_use, test_trial_to_use), [1, 2, 3]));
            n_BT_BI = length(control_trials_to_use);
        else
            BT_BI_mean = squeeze(mean(BT_BI_data.BT_ES(:, :, :, test_trial_to_use), [1, 2, 3]));
            n_BT_BI = length(control_trials_to_use_BT_BI);
        end

        n_BT_WI = length(control_trials_to_use);

        % Add new row to the summary table
        anat = unique(string(BT_WI_data.chan_test_table.anatomical_label));
        new_row = table(anat, WT_mean, BT_WI_mean, BT_BI_mean, n_BT_WI, n_BT_BI, test_trial_to_use, ...
                        'VariableNames', {'anat', 'mean_WT', 'mean_BT_WI', 'mean_BT_BI', 'n_BT_WI', 'n_BT_BI', 'test_trial_idx'});
        summary_table = [summary_table; new_row]; % Append the new row
    end
end

% Add difference column and save the table
summary_table.WIvsBI = summary_table.mean_BT_WI - summary_table.mean_BT_BI;

% Group summary for WIvsBI by anat
groupedby_anat_summary = groupsummary(summary_table, 'anat', {'mean'}, {'mean_WT', 'mean_BT_WI', 'mean_BT_BI', 'WIvsBI'});
groupedby_anat_summary = sortrows(groupedby_anat_summary, 'mean_WIvsBI', 'descend');

output_table_path = fullfile(params.output_destination, "WIvsBI_summary_table.mat");
save(output_table_path, "summary_table", "groupedby_anat_summary");

%% view average as a function of trial_idx to see if improved replay over time.
for chan_id = channel_ids_to_use

    % Load between_trials info for the channel (should be same whether
    % loaded from comparison{1} or comparison{2}

    % Load BT-BI data for this channel
    [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);
    %size(BT_ES) = 31   380    90    12 (enc, maint, ntrials_control, ntrials_test)

    flattenedMatrix = reshape(BT_BI_data.BT_ES, [], size(BT_BI_data.BT_ES, 3), size(BT_BI_data.BT_ES, 4));
    flattenedMatrix = squeeze(any(isnan(flattenedMatrix),1));
    test_trials_without_nan = ~all(flattenedMatrix, 1);
    ctrl_trials_without_nan = ~all(flattenedMatrix, 2);
    % figure;
    % imagesc(flattenedMatrix)

    % rows_without_nan = ~squeeze(any(isnan(BT_BI_data.BT_ES), [1,2,4]));
    BT_BI_data.BT_ES = BT_BI_data.BT_ES(:,:, ctrl_trials_without_nan, test_trials_without_nan);
    % rows_without_nan = ~squeeze(any(isnan(BT_WI_data.BT_ES), [1,2,4]));
    BT_WI_data.BT_ES = BT_WI_data.BT_ES(:,:, test_trials_without_nan, test_trials_without_nan); % here it is the same trials for dimensions 3,4
    BT_WI_data.chan_test_table = BT_WI_data.chan_test_table(test_trials_without_nan,:);

    % look at EMS as a function of between trials
    figure;
    num_trials = size(BT_WI_data.BT_ES, 4); % Number of lines to plot
    legend_entries = cell(1, num_trials);
    colors = jet(num_trials); % 'lines' colormap has distinct colors
    subplot(4,1,1) % trial_idx x axis
    for trial_idx = 1:num_trials
        plot(squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2])), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    ylim([0 0.8])
    xlabel('M trial_idx')
    subplot(4,1,2)
    for trial_idx = 1:num_trials
        plot(([1:num_trials] - trial_idx), squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2])), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    xlabel('E trial_idx - M trial_idx');
    ylim([0 0.8])
    subplot(4,1,4)
    for trial_idx = 1:num_trials
        plot(([BT_WI_data.chan_test_table.trial_ID(1:num_trials)] - BT_WI_data.chan_test_table.trial_ID(trial_idx)), squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2])), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    ylim([0 0.8])
    xlabel('E trial_id - M trial_id')
    subplot(4,1,3)
    for trial_idx = 1:num_trials
        plot(([BT_WI_data.chan_test_table.trial_ID(1:num_trials)]), squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2])), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    ylabel('mean EMS');
    xlabel('M trial_id')
    legend(legend_entries, 'Location', 'best');
    ylim([0 0.8])
    hold off;

        sgtitle([sprintf('EMS Mean Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
    patient_id, chan_id, enc_id, C{1, image_id}), ...
    unique(string(BT_WI_data.chan_test_table.anatomical_label))]);

    save_file = fullfile(params.output_destination, sprintf('BTWI_thru_trial_idx/Patient%d_Channel%d_Enc%d_Image%d.fig', ...
    patient_id, chan_id, enc_id, image_id));
    save_dir = fileparts(save_file);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    saveas(gcf, save_file);
    close()

        %% normalize each trace?
    figure;
    num_trials = size(BT_WI_data.BT_ES, 4); % Number of lines to plot
    legend_entries = cell(1, num_trials);
    colors = jet(num_trials); % 'lines' colormap has distinct colors
    subplot(4,1,1)
    for trial_idx = 1:num_trials
        plot(min_max_normalize(squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2]))), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    xlabel('trial_idx')
    subplot(4,1,2)
    for trial_idx = 1:num_trials
        plot(([1:num_trials] - trial_idx), min_max_normalize(squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2]))), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    xlabel('diff trial_idx');
    subplot(4,1,3)
    for trial_idx = 1:num_trials
        plot(([BT_WI_data.chan_test_table.trial_ID(1:num_trials)] - BT_WI_data.chan_test_table.trial_ID(trial_idx)), min_max_normalize(squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2]))), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    xlabel('diff trial_id')
    subplot(4,1,4)
    for trial_idx = 1:num_trials
        plot(([BT_WI_data.chan_test_table.trial_ID(1:num_trials)]), min_max_normalize(squeeze(mean(BT_WI_data.BT_ES(:,:,:,trial_idx), [1, 2]))), ...
             'Color', colors(trial_idx, :), 'LineWidth', 1.5); % Use unique color
        hold on;
        legend_entries{trial_idx} = ['Trial ', num2str(BT_WI_data.chan_test_table.trial_ID(trial_idx))];
    end
    xlabel('trial_id')
    ylabel('mean EMS');
    legend(legend_entries, 'Location', 'best');
    hold off;


    sgtitle([sprintf('Min-Max-Normalized Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
    patient_id, chan_id, enc_id, C{1, image_id}), ...
    unique(string(BT_WI_data.chan_test_table.anatomical_label))]);

    save_file = fullfile(params.output_destination, sprintf('BTWI_thru_trial_idx_MinMaxNormalized/Patient%d_Channel%d_Enc%d_Image%d.fig', ...
    patient_id, chan_id, enc_id, image_id));
    save_dir = fileparts(save_file);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    saveas(gcf, save_file);
    close()


end

function [y_values] = min_max_normalize(y_values)
    y_values = (y_values - min(y_values)) / (max(y_values) - min(y_values));
end

%% legacy
% %% average across pairs
%     % Subplot 2: Between-trial Within-Image similarity avg across pairs
%     subplot(1, 3, 2)
%     imagesc(squeeze(mean(BT_WI_data.BT_ES(:,:,:,test_trial_to_use),3,'omitnan')), [0 0.5])
%     title(strcat(comp_titles{2}, 'same trial vs WI avg across trial pairs n=', num2str(size(BT_WI_data.BT_ES(:,:,:,test_trial_to_use),3))))
% 
%     % Subplot 3: Between-trial Between-Image similarity avg across pairs
%     subplot(1, 3, 3)
%     imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,:, test_trial_to_use),3,'omitnan')), [0 0.5])
%     title(strcat(comp_titles{1}, 'same trial vs BI avg across trial pairs n=', num2str(size(BT_BI_data.BT_ES(:,:,:,test_trial_to_use),3))))
% 
% %% significance testing BT-WI and BT-BI
% %% TODO: put all together in 1 figure
% clearvars chan_test_table chan_id test_rows
% nPermutations=1000; % number of permutations for generation null distribution for cluster sizes
% use_z = false; % false uses p
% alpha = 0.05;
% z_thresh = 1.96;
% 
% addpath("../subfunctions")
% addpath("../clustering subfunctions")
% for chan_id = [10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66]
%     [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);
% 
%     figure;
%     hold on
%     subplot(2,2,1);
%     imagesc(squeeze(mean(BT_WI_data.BT_ES(:,:,:,:), [3 4], 'omitnan')))
%     clim([0 0.08])
%     title('BT WI mean across trial pairs and samples')
% 
%     hold on
%     subplot(2,2,2);
%     imagesc(squeeze(mean(BT_WI_data.BT_ES(:,:,1,1), [4], 'omitnan')))
%     clim([0 0.08])
%     title('BT WI single trial pair and sample')
% 
%     subplot(2,2,3);
%     imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,:,:),[3 4], 'omitnan')))
%     clim([0 0.08])
%     title('BT BI mean across trial pairs and samples')
% 
%     subplot(2,2,4);
%     imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,1,:),[4], 'omitnan')))
%     clim([0 0.08])
%     title('BT BI trial 1 pairs, meaned across samples')
% 
%     sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s, Anat: %s', ...
%         patient_id, chan_id, enc_id, C{1, image_id}, ...
%         unique(string(BT_WI_data.chan_test_table.anatomical_label)))]);
% 
%     savefig(fullfile(params.output_folder, ...
%         sprintf('Anat%s_Patient%d_Channel%d_Enc%d_Image%s.fig', ...
%         patient_id, chan_id, enc_id, C{1, image_id}, ...
%         unique(string(BT_WI_data.chan_test_table.anatomical_label)))))
% 
%     % Compute clusters
%     % for i =1:3 % BT_WI_data test sample id
%     %     [num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, ...
%     %      null_cluster_sizes, threshold_from_null, real_p_val_matrix] = ...
%     %      get_clusters(squeeze(mean(BT_BI_data,3)), ... % control mean x samples
%     %      squeeze(mean(BT_BI_data_std,3)), ... % control std x samples
%     %      squeeze(mean(BT_WI_data(:,:,:,i),[3])), ... % test mean % Enc,Maint, nTrialsWitem, 1000, samples
%     %      nPermutations, alpha, use_z, z_thresh);
%     % 
%     %     % visualize
%     % 
%     %     if max(max(real_p_val_matrix)) > 1-alpha
%     %         % error("test never signficant. greatest p:%s", num2str(max(max(real_p_val_matrix))))
%     %     else
%     %         fprintf("yay chan_id:%s\n", num2str(chan_id))
%     %         figure;
%     %         imagesc(real_p_val_matrix)
%     %         try
%     %             visualize_significant_clusters(real_p_val_matrix, sig_test_cluster_sizes, threshold_from_null, alpha); % view test sig/~sig clusters
%     %             sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
%     %                 patient_id, chan_id, enc_id, C{1, image_id}), ...
%     %                 unique(string(chan_test_table.anatomical_label))]);
%     %             plot_cluster_size_pdf(real_test_cluster_sizes, null_cluster_sizes, nPermutations); % test & null cluster size PDFs
%     %             sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
%     %                 patient_id, chan_id, enc_id, C{1, image_id}), ...
%     %                 unique(string(chan_test_table.anatomical_label))]);
%     %             xline(threshold_from_null, 'r', 'LineWidth', 2, 'DisplayName', sprintf('p < %d Threshold from Null', alpha)); % Significance testing
%     %         catch
%     %             fprintf("error on chan_id:%s\n",num2str(chan_id))
%     %         end
%     %     end
%     % end
% end
% 

% 
% %% ROI EMS vs EFS & 
% % recompute trial pairs
% % ROI Test EMS vs Ctrl EMS across all regions
% % recompute with original baselining and power computation.
% 
% %% check individual trial pairs & average across set of trial pairs
% figure;
% hold on
% subplot(1,3,1);
% imagesc(squeeze(mean(BT_WI_data.BT_ES(:,:,:,:), [3 4])))
% clim([0 0.08])
% subplot(1,3,2);
% imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,1,:),[4])))
% clim([0 0.08])
% subplot(1,3,3);
% imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,:,:),[3 4])))
% clim([0 0.08])
% 
% 
% %% checking wh ES_within is nans
% sum(isnan(all3_ES_matrix(:,1,1,1)))
% size(isnan(all3_ES_matrix(:,1,1,1)))
% 
% 
% sum(any((isnan(PSVs(:,:,:))|PSVs(:,:,:)==0), [1, 2]))
% size(isnan(PSVs(1,1,:)))
% 
% %%