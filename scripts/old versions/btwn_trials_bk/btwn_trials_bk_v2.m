%% v2 changing it so that BT-WI and BT-BI are subplots 2 and 3
%% v2 additional: adding sig. cluster analysis between BT-WI and BT-BI
% v1 WT-WI is subplot1, BT-WI or BT-BI is subplot 2.
%% between trials comparisons - within-item or between-items
% Set patient, image, and encoding IDs
patient_id = 201901;
image_id = 1;
enc_id = 2;
output_destination = 'C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBT & WIvsBI'; % Set your desired output folder

% Load required data
load(strcat("D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\",num2str(patient_id),"\\encoding_similarity.mat"))
load(strcat("D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\",num2str(patient_id),"\\PSVs.mat"), 'label_table')
load(strcat("D:\\Power Spectrum Similarity\\Z_Raw Data Storage\\OWM_trialinfo_",num2str(patient_id),".mat"), 'C')

rows_without_nan_by_enc = ~squeeze(any(any(isnan(all3_ES_matrix), 2), 3));

% Comparison options and titles
comp_options = {'corrWOitemMaint vs corrWitemEnc', 'corrWitemMaint vs corrWitemEnc'};
comp_titles = {'Between Trials - Between Images', 'Between Trials - Within Images'};

% Unique anatomical labels
[unique_channels_labels, ia] = unique(label_table.channel_ID, 'rows');
rows_without_nan = get_rows_without_nan(label_table);
label_table = label_table(rows_without_nan, :);

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
for chan_id = [10]%[10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66]

    % Load between_trials info for the channel (should be same whether
    % loaded from comparison{1} or comparison{2}
    load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\all_channels_data', ...
        num2str(patient_id), comp_options{2}, num2str(enc_id), num2str(image_id)))

    % Load BT-BI data for this channel
    load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\BT_%s.mat', ...
        num2str(patient_id), comp_options{1}, num2str(enc_id), num2str(image_id), num2str(chan_id)))
    BT_BI_data = control_EMS_matrices; % Rename for clarity
    BT_BI_data_std = control_EMS_matrices_std; % Rename for clarity

    % Load BT-WI data for this channel
    load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\BT_%s.mat', ...
        num2str(patient_id), comp_options{2}, num2str(enc_id), num2str(image_id), num2str(chan_id)))
    BT_WI_data = control_EMS_matrices; % Rename for clarity
    test_table = label_table(test_rows,:);
    chan_test_rows = test_table.channel_ID == chan_id; % Identify rows for this channel
    chan_test_table = test_table(chan_test_rows, :);

    Mtime = [250:640]; % Delay time

    figure;

    sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
        patient_id, chan_id, enc_id, C{1, image_id}), ...
        unique(string(chan_test_table.anatomical_label))]);

    % Subplot 1: Within-trial similarity
    subplot(1, 3, 1)
    result = all3_ES_matrix(rows_without_nan, :, Mtime, enc_id);
    result = result(test_rows, :, :, :);
    result = result(chan_test_rows, :, :, :);
    imagesc(squeeze(mean(result, 1, 'omitnan')), [0 0.5])
    title('Within-Trial Similarity')

    % Subplot 2: Between-trial Within-Image similarity
    subplot(1, 3, 2)
    imagesc(squeeze(mean(BT_WI_data, [3 4], 'omitnan')), [0 0.5])
    title(comp_titles{2})

    % Subplot 3: Between-trial Between-Image similarity
    subplot(1, 3, 3)
    imagesc(squeeze(mean(BT_BI_data, [3 4], 'omitnan')), [0 0.5])
    title(comp_titles{1})

    saveas(gcf, fullfile(output_destination, sprintf('Patient%d_Channel%d_Enc%d_Image%d.fig', ...
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
    [BT_BI_data, BT_BI_data_std, BT_WI_data, test_rows] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id);

    test_table = label_table(test_rows,:);
    chan_test_rows = test_table.channel_ID == chan_id; % Identify rows for this channel
    chan_test_table = test_table(chan_test_rows, :);

    figure;
    hold on
    subplot(2,2,1);
    imagesc(squeeze(mean(BT_WI_data(:,:,:,:), [4])))
    clim([0 0.08])
    title('BT WI mean across trial pairs and samples')
    subplot(2,2,4);
    imagesc(squeeze(mean(BT_BI_data(:,:,1,:),[4])))
    clim([0 0.08])
    title('BT BI trial 1 pairs, meaned across samples')
    subplot(2,2,3);
    imagesc(squeeze(mean(BT_BI_data(:,:,:,:),[3 4])))
    clim([0 0.08])
    title('BT BI mean across trial pairs and samples')

    sgtitle([sprintf('Patient: %d, Channel: %d, Enc: %d, Image: %s, Anat: %s', ...
        patient_id, chan_id, enc_id, C{1, image_id}, ...
        unique(string(chan_test_table.anatomical_label)))]);

    savefig(fullfile(output_folder, ...
        sprintf('Anat%s_Patient%d_Channel%d_Enc%d_Image%s', ...
        patient_id, chan_id, enc_id, C{1, image_id}, ...
        unique(string(chan_test_table.anatomical_label)))))

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

function [BT_BI_data, BT_BI_data_std, BT_WI_data, test_rows] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id)
    load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\all_channels_data', ...
        num2str(patient_id), comp_options{2}, num2str(enc_id), num2str(image_id)))

    % Load BT-BI data for this channel
    load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\BT_%s.mat', ...
        num2str(patient_id), comp_options{1}, num2str(enc_id), num2str(image_id), num2str(chan_id)))
    BT_BI_data = control_EMS_matrices; % Rename for clarity
    BT_BI_data_std = control_EMS_matrices_std; % Rename for clarity

    % Load BT-WI data for this channel
    load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\BT_%s.mat', ...
        num2str(patient_id), comp_options{2}, num2str(enc_id), num2str(image_id), num2str(chan_id)))
    BT_WI_data = control_EMS_matrices; % Rename for clarity
end

%% ROI EMS vs EFS & 
% recompute trial pairs
% ROI Test EMS vs Ctrl EMS across all regions
% recompute with original baselining and power computation.

%% check individual trial pairs & average across set of trial pairs
figure;
hold on
subplot(1,3,1);
imagesc(squeeze(mean(BT_WI_data(:,:,:,:), [4])))
clim([0 0.08])
subplot(1,3,2);
imagesc(squeeze(mean(BT_BI_data(:,:,1,:),[4])))
clim([0 0.08])
subplot(1,3,3);
imagesc(squeeze(mean(BT_BI_data(:,:,:,:),[3 4])))
clim([0 0.08])