%% EES item-by-item matrix
% for one patient, take all images enc1 and plot the average EES for that
% matrix vs all items
% v1 adapted from t_test_WI_vs_BI.m (or v1 of that script)

params = get_parameters();
params.btwn_trial_type = 'EES';
% Set patient, image, and encoding IDs
patient_id = 201901;
image_ids = 1:9;
enc_id = 1;
if patient_id == 201901
    channel_ids_to_use = [10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66];
elseif patient_id == 201908
    channel_ids_to_use = [28, 30, 37, 38, 74, 76, 82, 83, 84, 85];
elseif patient_id == 201910
    channel_ids_to_use = [7 14  25 26 27  41 42  52 53 54  59 60];
end

% Comparison options and titles
comp_options = {'corr BT BI', 'corr BT WI'};
comp_titles = {'Between Trials - Between Images', 'Between Trials - Within Images'};

% old=false; % whether to use true:(100 ms window, -.25--.75 baseline within trial) or false:(200 ms, -0.5-0.0 baseline across trials)
if strcmp(params.output_folder_name, 'allpatients gammamod allregions allitem allenc')
    params.output_destination = "C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBTvsBI newer_baselining_windowing"; % Set your desired output folder
elseif strcmp(params.output_folder_name, 'allpatients gammamod allregions allitem allenc baseline across trials')
    params.output_destination = 'C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBTvsBI original_baselining_windowing'; % Set your desired output folder
else
    ImplementationError("params.output_folder_name is not on the list")
end

% params.output_destination = fullfile(params.output_destination, sprintf("p%s_enc%s_image%s", num2str(patient_id), num2str(enc_id), num2str(image_id)));
params.output_destination = fullfile(params.output_destination, sprintf("p%s_enc%s_image%s", num2str(patient_id)));

% Load PSVs, within-trial similarity, images
load(fullfile(strcat(params.output_folder,"\\",num2str(patient_id),"\\encoding_similarity.mat"))) % within trial similarity
load(fullfile(strcat(params.output_folder,"\\",num2str(patient_id),"\\PSVs.mat")), 'label_table') % label table
patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C');
clear patient_preprocessed_data_path

function [BT_BI_data, BT_WI_data] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params)
    % Load BT-BI data for this channel
    BT_BI_data = load(sprintf('%s\\%s\\%s\\%s\\enc%s_image%s\\BT_%s.mat', ...
        params.output_folder, num2str(patient_id), comp_options{1}, params.btwn_trial_type, num2str(enc_id), num2str(image_id), num2str(chan_id)));
    clear BT_ES
    % BT_BI_data_std = control_EMS_matrices_std; % Rename for clarity

    % Load BT-WI data for this channel
    BT_WI_data = load(sprintf('%s\\%s\\%s\\%s\\enc%s_image%s\\BT_%s.mat', ...
        params.output_folder, num2str(patient_id), comp_options{2}, params.btwn_trial_type, num2str(enc_id), num2str(image_id), num2str(chan_id)));
end

function [logicalIdx1, logicalIdx2] = createLogicalIndex(sizeDim) % omits diagonal for WI pairs
    % createLogicalIndex creates two logical arrays for indexing a matrix,
    % excluding diagonal elements (e.g., 1,1; 2,2; ...; n,n).
    %
    % INPUT:
    %   sizeDim - Size of the dimension (e.g., 12 for a matrix (:,:,12,12))
    %
    % OUTPUT:
    %   logicalIdx1 - Logical index for the first dimension
    %   logicalIdx2 - Logical index for the second dimension

    % Validate input
    if ~isscalar(sizeDim) || sizeDim <= 0 || sizeDim ~= round(sizeDim)
        error('Input sizeDim must be a positive integer.');
    end

    % Initialize logical arrays
    logicalIdx1 = true(sizeDim, sizeDim);
    logicalIdx2 = true(sizeDim, sizeDim);

    % Create diagonal indices to exclude
    diagIdx = 1:sizeDim;

    % Set diagonal elements to false
    logicalIdx1(sub2ind([sizeDim, sizeDim], diagIdx, diagIdx)) = false;
    logicalIdx2(sub2ind([sizeDim, sizeDim], diagIdx, diagIdx)) = false;
end

function BT_struct = subset_by_all3_correct_test_trials(BT_struct)
    rows_with_all3_correct = sum(BT_struct.chan_test_table.encoding_correctness, 2)==3;
    BT_struct.chan_test_table = BT_struct.chan_test_table(rows_with_all3_correct,:);
    BT_struct.BT_ES = BT_struct.BT_ES(:,:,:,rows_with_all3_correct);
end

function BT_struct = subset_by_without_image_id(BT_struct, image_id)
    rows_without_image_id = all(BT_struct.chan_test_table.encID_to_imageID ~= image_id, 2);
    BT_struct.chan_test_table = BT_struct.chan_test_table(rows_without_image_id,:);
    BT_struct.BT_ES = BT_struct.BT_ES(:,:,:,rows_without_image_id);
end

%% main code

Ns_by_channel = nan(length(image_ids), length(image_ids), 2, length(channel_ids_to_use));
item_EES_matrix_by_channel = nan(length(image_ids), length(image_ids), length(channel_ids_to_use));
anat_by_channel = strings(1, length(channel_ids_to_use));
chan_idx=0; % iterator
for chan_id = channel_ids_to_use
    chan_idx=chan_idx+1;
    item_EES_matrix = nan(length(image_ids), length(image_ids));
    Ns = nan(length(image_ids), length(image_ids), 2);
    for first_image_id = image_ids
        for second_image_id = image_ids

            % get data
            [data_btwn_these_two_images, anat_this_channel] = get_data_btwn_these_two_images(first_image_id, second_image_id, patient_id, comp_options, enc_id, chan_id, params);

            % compute on this data
            Ns(first_image_id, second_image_id, :) = size(data_btwn_these_two_images, 3) * size(data_btwn_these_two_images, 4);
            mean_data_btwn_these_two_images = mean(data_btwn_these_two_images, [1:4], "omitnan"); % example bringing to single value
            item_EES_matrix(first_image_id, second_image_id) = mean_data_btwn_these_two_images;

            % store by channel
            Ns_by_channel(:,:,:,chan_idx) = Ns;
            item_EES_matrix_by_channel(:,:,chan_idx) = item_EES_matrix;
        end
    end
    anat_by_channel(chan_idx) = anat_this_channel;
end

%% plot figures
finite_values = item_EES_matrix_by_channel(~isinf(mean(item_EES_matrix_by_channel, 3, 'omitnan')));
mean_value = mean(finite_values); % Calculate mean and standard deviation
std_value = std(finite_values);
clims = [mean_value - 2 * std_value, mean_value + 2 * std_value]; % Set clims to 2 standard deviations (2z) away from the mean

figure;
title("average across channels")
imagesc(mean(item_EES_matrix_by_channel, 3, 'omitnan'))
set(gca, 'YTick', 1:length(C), 'YTickLabel', C); % Set the y-axis labels
set(gca, 'XTick', 1:length(C), 'XTickLabel', C); % Set the x-axis labels
xtickangle(45); % optional rotate x-axis labels to make them more visible
colorbar;
clim(clims)
sgtitle([sprintf('%s Patient: %d, Enc: %d', params.btwn_trial_type, patient_id, enc_id)]);



% clims = [min(item_EES_matrix_by_channel(:)), max(item_EES_matrix_by_channel(~isinf(item_EES_matrix_by_channel(:))))];
% clims = prctile(item_EES_matrix_by_channel(~isinf(item_EES_matrix_by_channel)), [5, 95]); % Calculate the 5th and 95th percentiles
finite_values = item_EES_matrix_by_channel(~isinf(item_EES_matrix_by_channel));
mean_value = mean(finite_values); % Calculate mean and standard deviation
std_value = std(finite_values);
clims = [mean_value - 1 * std_value, mean_value + 1 * std_value]; % Set clims to 2 standard deviations (2z) away from the mean

figure;
num_channels = length(channel_ids_to_use);
len_dim = ceil(sqrt(num_channels));
for chan_idx = 1:num_channels
    subplot(len_dim, len_dim, chan_idx)
    imagesc(item_EES_matrix_by_channel(:,:,chan_idx))
    title(strcat(anat_by_channel(chan_idx), num2str(channel_ids_to_use(chan_idx))))
    colorbar;
    clim(clims)
end
sgtitle([sprintf('%s Patient: %d, Enc: %d', params.btwn_trial_type, patient_id, enc_id)]);
% xlabel(gcf, 'Image ID', 'FontSize', 14);
% ylabel(gcf, 'Image ID', 'FontSize', 14);
% ylabel(colorbar, 'rho');%, 'FontSize', 12, 'Rotation', 270, 'VerticalAlignment', 'bottom');
%% plotting N
Ns_by_channel_reduced = squeeze(prod(Ns_by_channel, 3));

mean_value = mean(sum(Ns_by_channel_reduced, 3), [1:4]); % Calculate mean and standard deviation
std_value = mean(std(sum(Ns_by_channel_reduced, 3)));
clims = [mean_value - 2 * std_value, mean_value + 2 * std_value]; % Set clims to 2 standard deviations (2z) away from the mean

figure;
title("N pairs for corresponding plot")
imagesc(sum(Ns_by_channel_reduced, 3, 'omitnan'))
set(gca, 'YTick', 1:length(C), 'YTickLabel', C); % Set the y-axis labels
set(gca, 'XTick', 1:length(C), 'XTickLabel', C); % Set the x-axis labels
xtickangle(45); % optional rotate x-axis labels to make them more visible
colorbar;
% clim(clims)
sgtitle([sprintf('Ns %s Patient: %d, Enc: %d', params.btwn_trial_type, patient_id, enc_id)]);


figure;
num_channels = length(channel_ids_to_use);
len_dim = ceil(sqrt(num_channels));
for chan_idx = 1:num_channels
    subplot(len_dim, len_dim, chan_idx)
    imagesc(Ns_by_channel_reduced(:,:,chan_idx))
    title(strcat(anat_by_channel(chan_idx), num2str(channel_ids_to_use(chan_idx))))
    colorbar;
    % clim(clims)
end
sgtitle([sprintf('Ns %s Patient: %d, Enc: %d', params.btwn_trial_type, patient_id, enc_id)]);
%% troubleshooting
chan_id = 22;
patient_id = 201901;
[first_image_id, second_image_id] = deal(3,4);
data_btwn_these_two_images_1 = get_data_btwn_these_two_images(first_image_id, second_image_id, patient_id, comp_options, enc_id, chan_id, params);
[first_image_id, second_image_id] = deal(4,3);
data_btwn_these_two_images_2 = get_data_btwn_these_two_images(first_image_id, second_image_id, patient_id, comp_options, enc_id, chan_id, params);

disp(sum(data_btwn_these_two_images_1 == data_btwn_these_two_images_2), [1:4]);
disp(size(data_btwn_these_two_images_1))


% % Compute the flattened size of the first two dimensions
% size_flattened = size(data_btwn_these_two_images_1, 1) * size(data_btwn_these_two_images_1, 2);
% 
% % Reshape the matrix
% flattened_matrix = reshape(data_btwn_these_two_images_1 == data_btwn_these_two_images_2_permuted, ...
%                            size_flattened, size(data_btwn_these_two_images_1, 3), size(data_btwn_these_two_images_1, 4));

function [data_btwn_these_two_images, anat] = get_data_btwn_these_two_images(first_image_id, second_image_id, patient_id, comp_options, enc_id, chan_id, params)
    % get data
    [BI.EES, WI.EES] = load_BT_data(patient_id, comp_options, enc_id, first_image_id, chan_id, params);
    BI.EES = subset_by_all3_correct_test_trials(BI.EES); % remove trials that are not all3 correct to make the data symmetrical
    BI.EES = subset_by_without_image_id(BI.EES, second_image_id);  % remove trials where the second_image_id will appear to make the data symmetrical
    rows_with_second_image = BI.EES.chan_ctrl_table.encID_to_imageID(:,enc_id) == second_image_id;
    if first_image_id == second_image_id
        if sum(rows_with_second_image) ~= 0
            error("'between image' data indicates that some 'within image' data is present")
        end
        [~, logicalIdx2] = createLogicalIndex(size(WI.EES.BT_ES, 3)); % indices that omit same-trial pairs for WI.
        data_btwn_these_two_images = WI.EES.BT_ES(:,:,logicalIdx2);
    else
        data_btwn_these_two_images = BI.EES.BT_ES(:,:,rows_with_second_image,:);
        data_btwn_these_two_images = data_btwn_these_two_images(:,1:31,:,:);    % convert from all3encoding windows to the first one
    end

    anat = unique(string(BI.EES.chan_ctrl_table.anatomical_label));
end