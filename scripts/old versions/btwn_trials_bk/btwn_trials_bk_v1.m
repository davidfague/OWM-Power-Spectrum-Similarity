%% v2 changing it so that BT-WI and BT-BI are subplots 2 and 3
% v1 WT-WI is subplot1, BT-WI or BT-BI is subplot 2.
%% between trials comparions - within-item or between-items

load("D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem allenc\201907\encoding_similarity.mat")
load("D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem allenc\201907\PSVs.mat", 'label_table')
% load("D:\Power Spectrum Similarity\Z_Raw Data Storage\D_OWM_t_bipolar_201907.mat")
load("D:\Power Spectrum Similarity\Z_Raw Data Storage\OWM_trialinfo_201907.mat", 'C')
rows_without_nan_by_enc = ~squeeze(any(any(isnan(all3_ES_matrix), 2), 3));

%%
patient_id = 201910;
image_id = 1;
enc_id = 1;

comp_options = {'corrWOitemMaint vs corrWitemEnc', 'corrWitemMaint vs corrWitemEnc'};
comparison = 1; % 1 for btwn-image ; 2 for within-image
if comparison==1
    comp_title = 'Between Images';
    comparison = comp_options(1);
elseif comparison==2
    comp_title = 'Within Images';
    comparison = comp_options(2);
else
    error("comparison not implemented")
end

% load images
load(strcat("D:\Power Spectrum Similarity\Z_Raw Data Storage\OWM_trialinfo_",num2str(patient_id),".mat"), 'C')

% load within trial data
load(strcat("D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem allenc\",num2str(patient_id),"\encoding_similarity.mat"))

% load between_trials info
load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\all_channels_data', ...
    num2str(patient_id), comparison{1}, num2str(enc_id), num2str(image_id)))

% load label_table
load(strcat("D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem allenc\", num2str(patient_id), "\PSVs.mat"), 'label_table')
[unique_channels_labels, ia] = unique(label_table.channel_ID, 'rows');
% disp(label_table(ia,2:3));

rows_without_nan = get_rows_without_nan(label_table);
test_table = label_table(rows_without_nan,:);
test_table = test_table(test_rows,:);
%%
% amyg, hipp channels
channels_to_use = contains(string(label_table{ia, 3}), 'amyg', 'IgnoreCase', true) | ...
                  contains(string(label_table{ia, 3}), 'hipp', 'IgnoreCase', true);
channels_to_use = label_table{ia(channels_to_use),2};
% for chan_id = [10 11 12  22  61] % 201901
% for chan_id = [1 2  10 11 12 13  46 47 48 49 50  57 58 66] % 201910
%% download amyg, hipp BT data
%% download BI data

% mfg, mtg channels
% for chan_id = [7 8 9   22 23 24 25 26   48 49 50   63 64 65 66] % 201901
% for chan_id = [7 14  25 26 27  41 42  52 53 54  59 60] % 201910
% for chan_id = [90 18 26 28 36 39] % was for 201907 but weird data
% chan_id = 39;
for chan_id = [10 11 12  22  61 7 8 9   22 23 24 25 26   48 49 50   63 64 65 66]
% load between_trials data
load(sprintf('D:\\Power Spectrum Similarity\\AA_Processed Data\\allpatients gammamod allregions allitem allenc\\%s\\%s\\between_trials_enc%s_image%s\\BT_%s.mat', ...
    num2str(patient_id), comparison{1}, num2str(enc_id), num2str(image_id), num2str(chan_id)))
%%

chan_test_rows = test_table.channel_ID == chan_id; % Identify rows for this channel
chan_test_table = test_table(chan_test_rows,:);

Mtime=[250:640]; % delay time

figure;

sgtitle([comp_title, ...
         patient_id, ...
         strcat('enc', num2str(enc_id)), ...
         C{1, image_id}, ...
         num2str(chan_id), ...
         unique(string(chan_test_table.anatomical_label))]);

% within trial
subplot(1,2,1)
result = all3_ES_matrix(rows_without_nan, :, Mtime, enc_id);
result = result(test_rows, :, :, :);
result = result(chan_test_rows, :, :, :);
imagesc(squeeze(mean(result, 1, 'omitnan')), [0 0.5])
title('within trial similarity')

% btwn trial 
subplot(1,2,2)
imagesc(squeeze(mean(control_EMS_matrices, [3 4], 'omitnan')), [0 0.5])
title('between trial similarity')
end