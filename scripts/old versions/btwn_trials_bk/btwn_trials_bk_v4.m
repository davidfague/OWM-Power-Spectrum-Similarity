%% between trials comparisons - within-item and between-items
% v4 update file names,
% v3 enabling params; remove nan filtering based on table.
% v2 changing it so that BT-WI and BT-BI are subplots 2 and 3
% v2 additional: adding sig. cluster analysis between BT-WI and BT-BI
% v1 WT-WI is subplot1, BT-WI or BT-BI is subplot 2.

params = get_parameters();

% Set patient, image, and encoding IDs
patient_id = 201901;
image_id = 1;
enc_id = 1;

% Comparison options and titles
comp_options = {'corr BT BI', 'corr BT WI'};
comp_titles = {'Between Trials - Between Images', 'Between Trials - Within Images'};

old=true; % whether to use true:(100 ms window, -.25--.75 baseline within trial) or false:(200 ms, -0.5-0.0 baseline across trials)
if old
    output_folder_name = 'allpatients gammamod allregions allitem allenc';
    params.output_destination = "C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBT & WIvsBI ('newer' baselining,windowing) (still not all3enc set)"; % Set your desired output folder
else
    output_folder_name = 'allpatients gammamod allregions allitem allenc baseline across trials';
    params.output_destination = 'C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBT & WIvsBI (original baselining,windowing) (still not all3enc set)'; % Set your desired output folder
end
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

%% plot
[all_channels_data] = load_all_channels_data(params, patient_id, comp_options, enc_id, image_id);

for chan_id = [10]%[10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66]

    % Load between_trials info for the channel (should be same whether
    % loaded from comparison{1} or comparison{2}

    % Load BT-BI data for this channel
    [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);

    %size(BT_ES) = 31   380    90    12 (enc, maint, ntrials_test, ntrials_control)


    Mtime = [250:630]; % Delay time

    figure;

    sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
        patient_id, chan_id, enc_id, C{1, image_id}), ...
        unique(string(BT_WI_data.chan_test_table.anatomical_label))]);

    % Subplot 1: Within-trial similarity for trial 1 of test
    subplot(1, 3, 1)
    test_trial_to_use = 1; % first, second, etc.; change if nan % test refers to target_item_correct_trials
    original_row_for_test_trial_1 = find((label_table.trial_ID == BT_BI_data.chan_test_table(test_trial_to_use,:).trial_ID) & (label_table.channel_ID == BT_BI_data.chan_test_table(test_trial_to_use,:).channel_ID));
    within_trial_single_trial = squeeze(all3_ES_matrix(original_row_for_test_trial_1, :, Mtime, enc_id)); % result = all3_ES_matrix(rows_without_nan, :, Mtime, enc_id);
    if any(isnan(within_trial_single_trial(:)))
        error("nan in this within trial data. change test_trial_to_use")
    end
    % result = result(test_rows, :, :, :);
    % result = result(chan_test_rows, :, :, :);
    % imagesc(squeeze(mean(result, 1, 'omitnan')), [0 0.5])
    imagesc(within_trial_single_trial, [0 0.5])
    title('Within-Trial Similarity Mean across trial')

    % Subplot 2: Between-trial Within-Image similarity
    subplot(1, 3, 2)
    imagesc(squeeze(mean(BT_WI_data.BT_ES, [3 4], 'omitnan')), [0 0.5])
    title(strcat(comp_titles{2}, ' mean across trial pairs'))

    % Subplot 3: Between-trial Between-Image similarity
    subplot(1, 3, 3)
    imagesc(squeeze(mean(BT_BI_data.BT_ES, [3 4], 'omitnan')), [0 0.5])
    title(comp_titles{1}, ' mean across trial pairs')

    if ~exist(params.output_destination, 'dir')
        mkdir(params.output_destination)
    end
    saveas(gcf, fullfile(params.output_destination, sprintf('Patient%d_Channel%d_Enc%d_Image%d.fig', ...
    patient_id, chan_id, enc_id, image_id)));

end

%% significance testing BT-WI and BT-BI
%% TODO: put all together in 1 figure
clearvars chan_test_table chan_id test_rows
nPermutations=1000; % number of permutations for generation null distribution for cluster sizes
use_z = false; % false uses p
alpha = 0.05;
z_thresh = 1.96;

addpath("../subfunctions")
addpath("../clustering subfunctions")
for chan_id = [10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66]
    [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);

    figure;
    hold on
    subplot(2,2,1);
    imagesc(squeeze(mean(BT_WI_data.BT_ES(:,:,:,:), [3 4], 'omitnan')))
    clim([0 0.08])
    title('BT WI mean across trial pairs and samples')

    hold on
    subplot(2,2,2);
    imagesc(squeeze(mean(BT_WI_data.BT_ES(:,:,1,1), [4], 'omitnan')))
    clim([0 0.08])
    title('BT WI single trial pair and sample')

    subplot(2,2,3);
    imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,:,:),[3 4], 'omitnan')))
    clim([0 0.08])
    title('BT BI mean across trial pairs and samples')

    subplot(2,2,4);
    imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,1,:),[4], 'omitnan')))
    clim([0 0.08])
    title('BT BI trial 1 pairs, meaned across samples')

    sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s, Anat: %s', ...
        patient_id, chan_id, enc_id, C{1, image_id}, ...
        unique(string(BT_WI_data.chan_test_table.anatomical_label)))]);

    savefig(fullfile(params.output_folder, ...
        sprintf('Anat%s_Patient%d_Channel%d_Enc%d_Image%s.fig', ...
        patient_id, chan_id, enc_id, C{1, image_id}, ...
        unique(string(BT_WI_data.chan_test_table.anatomical_label)))))

    % Compute clusters
    % for i =1:3 % BT_WI_data test sample id
    %     [num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, ...
    %      null_cluster_sizes, threshold_from_null, real_p_val_matrix] = ...
    %      get_clusters(squeeze(mean(BT_BI_data,3)), ... % control mean x samples
    %      squeeze(mean(BT_BI_data_std,3)), ... % control std x samples
    %      squeeze(mean(BT_WI_data(:,:,:,i),[3])), ... % test mean % Enc,Maint, nTrialsWitem, 1000, samples
    %      nPermutations, alpha, use_z, z_thresh);
    % 
    %     % visualize
    % 
    %     if max(max(real_p_val_matrix)) > 1-alpha
    %         % error("test never signficant. greatest p:%s", num2str(max(max(real_p_val_matrix))))
    %     else
    %         fprintf("yay chan_id:%s\n", num2str(chan_id))
    %         figure;
    %         imagesc(real_p_val_matrix)
    %         try
    %             visualize_significant_clusters(real_p_val_matrix, sig_test_cluster_sizes, threshold_from_null, alpha); % view test sig/~sig clusters
    %             sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
    %                 patient_id, chan_id, enc_id, C{1, image_id}), ...
    %                 unique(string(chan_test_table.anatomical_label))]);
    %             plot_cluster_size_pdf(real_test_cluster_sizes, null_cluster_sizes, nPermutations); % test & null cluster size PDFs
    %             sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
    %                 patient_id, chan_id, enc_id, C{1, image_id}), ...
    %                 unique(string(chan_test_table.anatomical_label))]);
    %             xline(threshold_from_null, 'r', 'LineWidth', 2, 'DisplayName', sprintf('p < %d Threshold from Null', alpha)); % Significance testing
    %         catch
    %             fprintf("error on chan_id:%s\n",num2str(chan_id))
    %         end
    %     end
    % end
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

%% ROI EMS vs EFS & 
% recompute trial pairs
% ROI Test EMS vs Ctrl EMS across all regions
% recompute with original baselining and power computation.

%% check individual trial pairs & average across set of trial pairs
figure;
hold on
subplot(1,3,1);
imagesc(squeeze(mean(BT_WI_data.BT_ES(:,:,:,:), [3 4])))
clim([0 0.08])
subplot(1,3,2);
imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,1,:),[4])))
clim([0 0.08])
subplot(1,3,3);
imagesc(squeeze(mean(BT_BI_data.BT_ES(:,:,:,:),[3 4])))
clim([0 0.08])


%% checking wh ES_within is nans
sum(isnan(all3_ES_matrix(:,1,1,1)))
size(isnan(all3_ES_matrix(:,1,1,1)))


sum(any((isnan(PSVs(:,:,:))|PSVs(:,:,:)==0), [1, 2]))
size(isnan(PSVs(1,1,:)))

%%