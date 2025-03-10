%%TODO:
% add N x N,
% add t-mask all-patient by-region analysis.
%%
% plot raw & t-test between WI and BI

% later add WC vs BC

% matching Fig. 2D of Pacheco-Estafan et al 2023

% notes: They used windows with (size, step)=(500,10) we use (200,10)

% we will be looking at frequencies 1:40 instead of bands 3-8, 9-12, 13-29,
% 3-75, 75-150.

% they include all electrodes in the region into the representational
% feature vector while we include only 1 electrode at a time.

% We look at gamma modulated channels specifically (can calculate non-gamma
% as well)

% it appears they fisher z-transform their rho values.

% loading can be done similar to btwn_trial_bk_v6.m

% in their code:
% https://github.com/dpachec/WM/blob/main/final_analysis_figures_3.m
% they use ttest() to ttest and atanh() to Fisher Z transform.

params = get_parameters();
% params.btwn_trial_type = 'EES';
% Set patient, image, and encoding IDs
% patient_id = 201908;
image_id = 1;
enc_id = 1;
% if patient_id == 201901
%     channel_ids_to_use = [10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66];
% elseif patient_id == 201908
%     channel_ids_to_use = [28, 30, 37, 38, 74, 76, 82, 83, 84, 85];
% elseif patient_id == 201910
%     channel_ids_to_use = [7 14  25 26 27  41 42  52 53 54  59 60];
% end

addpath("C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\subfunctions")
% Comparison options and titles
comp_options = {'corr BT BI', 'corr BT WI'};
comp_titles = {'Between Trials - Between Images', 'Between Trials - Within Images'};

% old=false; % whether to use true:(100 ms window, -.25--.75 baseline within trial) or false:(200 ms, -0.5-0.0 baseline across trials)
% if old
%     output_folder_name = 'allpatients gammamod allregions allitem allenc';
%     params.output_destination = "C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBTvsBI newer_baselining_windowing"; % Set your desired output folder
% else
%     output_folder_name = 'allpatients gammamod allregions allitem allenc baseline across trials';
%     params.output_destination = 'C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\Results\WTvsBTvsBI original_baselining_windowing'; % Set your desired output folder
% end

% params.output_destination = fullfile(params.output_destination, sprintf("p%s_enc%s_image%s", num2str(patient_id), num2str(enc_id), num2str(image_id)));
% clear old

% Load PSVs, within-trial similarity, images
% load(fullfile(strcat(params.output_folder,"\\",num2str(patient_id),"\\encoding_similarity.mat"))) % within trial similarity
% load(fullfile(strcat(params.output_folder,"\\",num2str(patient_id),"\\PSVs.mat")), 'label_table') % label table
% patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
% load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C');
% clear patient_preprocessed_data_path

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

function save_p_values(patient_id, comp_options, enc_id, image_id, params, p_values_by_chan, anat_by_chan, chan_id_by_chan)
    save(sprintf('%s\\%s\\%s\\%s\\enc%s_image%s\\p_values.mat',  params.output_folder, num2str(patient_id), comp_options{2}, params.btwn_trial_type, num2str(enc_id), num2str(image_id)), ...
        'p_values_by_chan', 'anat_by_chan', 'chan_id_by_chan')
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

%% plot raw data single & average EMS, EES across WI pairs and BI pairs
function plot_row(last_subplot_id, trial_idx1, trial_idx2, rho_lims, WI, BI)
    subplot(4,4,last_subplot_id+1) %EES
    imagesc(WI.EES.BT_ES(:,:,trial_idx1,trial_idx2))
    title('WI EES single pair')
    clim(rho_lims)

    subplot(4,4,last_subplot_id+2)
    imagesc(WI.EMS.BT_ES(:,:,trial_idx1,trial_idx2))
    title('WI EMS single pair')
    clim(rho_lims)

    subplot(4,4,last_subplot_id+3)
    imagesc(BI.EES.BT_ES(:,:,trial_idx1,trial_idx2))
    title('BI EES single pair')
    xline([40, 90, 130], 'LineStyle', '--')
    clim(rho_lims)

    subplot(4,4,last_subplot_id+4)
    imagesc(BI.EMS.BT_ES(:,:,trial_idx1,trial_idx2))
    title(sprintf('BI EMS ctrl=%d, test=%d', trial_idx1, trial_idx2))
    clim(rho_lims)
end

for chan_id = channel_ids_to_use(1)

    % Load between_trials info for the channel (should be same whether
    % loaded from comparison{1} or comparison{2}

    % Load BT-BI data for this channel
    WI = struct();
    BI = struct();
    params.btwn_trial_type = 'EES';
    [BI.EES, WI.EES] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);
    params.btwn_trial_type = 'EMS';
    [BI.EMS, WI.EMS] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);

    [logicalIdx1, logicalIdx2] = createLogicalIndex(size(WI.EMS.BT_ES, 3)); % indices that omit same-trial pairs for WI.

    if strcmp(params.btwn_trial_type, 'EMS')
        Mtime = [250:630]; % Delay time
    elseif strcmp(params.btwn_trial_type, 'EES')
        Mtime = [100:150];
    end

    figure;
    sgtitle([sprintf('%s Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
        params.btwn_trial_type, patient_id, chan_id, enc_id, C{1, image_id}), ...
        unique(string(WI.EMS.chan_test_table.anatomical_label))]);

    rho_lims = [0 0.05];

    % row1: averages;      WI:EES(1),EMS(2) BI:EES(3),EMS(4) %1-4
    subplot(4,4,1) %EES
    imagesc(mean(WI.EES.BT_ES(:,:,logicalIdx2), 3, 'omitnan'))
    title('WI EES avg')
    clim(rho_lims)

    subplot(4,4,2)
    imagesc(mean(WI.EMS.BT_ES(:,:,logicalIdx2), 3, 'omitnan'))
    title('WI EMS avg')
    clim(rho_lims)

    subplot(4,4,3)
    imagesc(mean(BI.EES.BT_ES(:,:,:,:), [3 4], 'omitnan'))
    title('BI EES avg')
    xline([40, 90, 130], 'LineStyle', '--')
    clim(rho_lims)

    subplot(4,4,4)
    imagesc(mean(BI.EMS.BT_ES(:,:,:,:), [3 4], 'omitnan'))
    title('BI EMS avg')
    clim(rho_lims)

    % other rows: individual trial pairs
    last_row = 4;
    for trial_idx1=1:size(WI.EMS.BT_ES, 3)
        for trial_idx2=1:size(WI.EMS.BT_ES, 3)
            if trial_idx1 == trial_idx2 || any(isnan(WI.EMS.BT_ES(:,:,trial_idx1, trial_idx2)), [1 2]) || any(isnan(BI.EMS.BT_ES(:,:,trial_idx1, trial_idx2)), [1 2]) || last_row == 16
                continue % skip same pairs & nans
            end
            

            plot_row(last_row, trial_idx1, trial_idx2, rho_lims, WI, BI) % row2: trial1:trial2; WI:EES(1),EMS(2) BI:EES(3),EMS(4) %5-8
            last_row = last_row +4;
            % plot_row(8, 2, 3, rho_lims, WI, BI) % row3: trial2:trial3; WI:EES(1),EMS(2) BI:EES(3),EMS(4) %9-12
            % 
            % plot_row(12, 1, 3, rho_lims, WI, BI) % row4: trial1:trial3; WI:EES(1),EMS(2) BI:EES(3),EMS(4) % 13-16
        end
    end
end
%% plot t-test between WI and BI.
color_limits = [-7 7];
function diff = calc_diff_all_pairs(WI, BI, num_pairs)
    % Inputs:
    % WI - The within-item matrix (size: time x time x num_WI)
    % BI - The between-item matrix (size: time x time x num_BI)
    % num_pairs - Maximum number of pairs to compute (randomly selected if less than total_pairs)

    % Dimensions of the input matrices
    num_WI = size(WI, 3);
    num_BI = size(BI, 3);

    % Determine the total number of pairs
    total_pairs = num_WI * num_BI;

    % Validate num_pairs
    if nargin < 3 || isempty(num_pairs)
        num_pairs = total_pairs; % Default to all pairs if num_pairs not specified
    end
    num_pairs = min(num_pairs, total_pairs);

    % Generate random pair indices if num_pairs is less than total_pairs
    if num_pairs < total_pairs
        % rng(0); % For reproducibility (optional, remove if not needed)
        random_indices = randperm(total_pairs, num_pairs);
    else
        random_indices = 1:total_pairs; % Use all combinations if num_pairs equals total_pairs
    end

    % Preallocate the output difference matrix
    diff = zeros(size(WI, 1), size(WI, 2), num_pairs);

    % Calculate differences for the selected combinations
    for k = 1:num_pairs
        [i, j] = ind2sub([num_WI, num_BI], random_indices(k)); % Map linear index to subscripts
        diff(:, :, k) = WI(:, :, i) - BI(:, :, j);
    end
end


function [t_values, p_values] = calc_ttest(diff)
    dim1 = size(diff,1);
    dim2 = size(diff,2);
    t_values = zeros(dim1, dim2);
    p_values = zeros(dim1, dim2);
    for i = 1:dim1
        for j = 1:dim2
            % Extract the distribution for the current time-time pair
            data_diff = squeeze(diff(i, j, :));
            
            % Perform one-sample t-test against 0
            [~, p, ~, stats] = ttest(data_diff, 0);
            
            % Store the t-value and p-value
            t_values(i, j) = stats.tstat;
            p_values(i, j) = p;
        end
    end
end

for chan_id = channel_ids_to_use

    % Load between_trials info for the channel (should be same whether
    % loaded from comparison{1} or comparison{2}

    % Load BT-BI data for this channel
    WI = struct();
    BI = struct();
    params.btwn_trial_type = 'EES';
    [BI.EES, WI.EES] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);
    params.btwn_trial_type = 'EMS';
    [BI.EMS, WI.EMS] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);

    [logicalIdx1, logicalIdx2] = createLogicalIndex(size(WI.EMS.BT_ES, 3)); % indices that omit same-trial pairs for WI.

    if strcmp(params.btwn_trial_type, 'EMS')
        Mtime = [250:630]; % Delay time
    elseif strcmp(params.btwn_trial_type, 'EES')
        Mtime = [100:150];
    end
    figure;
    sgtitle([sprintf('%s Patient: %d, Channel: %d, Enc: %d, Image: %s', ...
        params.btwn_trial_type, patient_id, chan_id, enc_id, C{1, image_id}), ...
        unique(string(WI.EMS.chan_test_table.anatomical_label))]);
    
    % EES_WI = mean(WI.EES.BT_ES(:,:,logicalIdx2), 3, 'omitnan');
    % EES_BI = mean(BI.EES.BT_ES(:,[1:31],:,:), [3 4], 'omitnan');
    % EES_diff = EES_WI - EES_BI;
    EES_WI = WI.EES.BT_ES(:,:,logicalIdx2);
    EES_BI = BI.EES.BT_ES(:,[1:31],:);
    EES_BI  = EES_BI(:,:,all(~isnan(EES_BI),[1 2]));%remove nans trial pairs
    EES_WI  = EES_WI(:,:,all(~isnan(EES_WI),[1 2]));%remove nans trial pairs
    
    EMS_WI = WI.EMS.BT_ES(:,:,logicalIdx2);
    EMS_BI = BI.EMS.BT_ES(:,:,:);
    EMS_WI  = EMS_WI(:,:,all(~isnan(EMS_WI),[1 2]));%remove nans trial pairs
    EMS_BI  = EMS_BI(:,:,all(~isnan(EMS_BI),[1 2]));%remove nans trial pairs
    
    
    EES_diff = calc_diff_all_pairs(EES_WI, EES_BI, 10000);
    EMS_diff = calc_diff_all_pairs(EMS_WI, EMS_BI, 10000);
    
    [EES_t_values, EES_p_values] = calc_ttest(EES_diff);
    [EMS_t_values, EMS_p_values] = calc_ttest(EMS_diff);
    
    subplot(2,4,1)
    imagesc(mean(EES_WI, 3, 'omitnan'))
    title('EES WI avg')
    colorbar;
    
    subplot(2,4,2)
    imagesc(mean(EES_BI, 3, 'omitnan'))
    title('EES BI avg')
    colorbar;
    
    subplot(2,4,3)
    imagesc(mean(EES_diff, 3, 'omitnan'))
    title('EES WI-BI avg')
    colorbar;
    
    subplot(2,4,4)
    imagesc(EES_t_values)
    title('EES t(WI-BI,0)')
    clim(color_limits);
    colormap('hot');
    colormap(flipud(colormap));
    colorbar;
    
    subplot(2,4,5)
    imagesc(mean(EMS_WI, 3, 'omitnan'))
    title('EMS WI avg')
    colorbar;
    
    subplot(2,4,6)
    imagesc(mean(EMS_BI, 3, 'omitnan'))
    title('EMS BI avg')
    colorbar;
    
    subplot(2,4,7)
    imagesc(mean(EMS_diff, 3, 'omitnan'))
    title('EMS WI-BI avg')
    colorbar;
    
    subplot(2,4,8)
    imagesc(EMS_t_values)
    title('EMS t(WI-BI,0)')
    colormap('hot')
    colormap(flipud(colormap));
    clim(color_limits)
    colorbar;

    if mean(EES_t_values, [1 2]) > 7
        tag = '_pos';
    elseif mean(EES_t_values, [1 2]) < -7
        tag = '_neg';
    else
        tag = '';
    end

    save_file = fullfile(params.output_destination, sprintf('BTWI_t_testing/t_score_visualization/Patient%d_Channel%d_Enc%d_Image%d%s.fig', ...
    patient_id, chan_id, enc_id, image_id, tag));
    save_dir = fileparts(save_file);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    saveas(gcf, save_file);
end

%% do significance testing using 'surrogate'
for patient_id = params.patient_IDs
   patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
   load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C');
   load(fullfile(patient_preprocessed_data_path, "D_OWM_t.mat"), 'labelsanatbkedit');
   channel_ids_to_use = load_channels_from_file(patient_id, comp_options, enc_id, image_id, params);
   original_labels = labelsanatbkedit.anatmacro1(channel_ids_to_use);
   channel_ids_to_use = filter_channels_ids_to_use(channel_ids_to_use, original_labels);

   anat_by_channels = labelsanatbkedit.anatmacro1(channel_ids_to_use);
    for image_id = 1:9
        final_p_by_channels_image = nan(length(channel_ids_to_use), 9);
        % anat_by_channels = strings(length(channel_ids_to_use), 1);
    
        if exist(sprintf("p_values_by_channel_%s.mat", num2str(patient_id)), 'file')
            load(sprintf("p_values_by_channel_%s.mat", num2str(patient_id)))
        else
            % add image_id dimension to final_ps
            for chan_idx = 1:length(channel_ids_to_use)
                chan_id = channel_ids_to_use(chan_idx);
            
                % Load BT-BI data for this channel
                WI = struct();
                BI = struct();
                params.btwn_trial_type = 'EES';
                [BI.EES, WI.EES] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);
                params.btwn_trial_type = 'EMS';
                [BI.EMS, WI.EMS] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params);
            
                [logicalIdx1, logicalIdx2] = createLogicalIndex(size(WI.EMS.BT_ES, 3)); % indices that omit same-trial pairs for WI.
            
                if strcmp(params.btwn_trial_type, 'EMS')
                    Mtime = [250:630]; % Delay time
                elseif strcmp(params.btwn_trial_type, 'EES')
                    Mtime = [100:150];
                end
            
                % EES_diff = EES_WI - EES_BI;
                EES_WI = WI.EES.BT_ES(:,:,logicalIdx2);
                EES_BI = BI.EES.BT_ES(:,[1:31],:);
                EES_BI  = EES_BI(:,:,all(~isnan(EES_BI),[1 2]));%remove nans trial pairs
                EES_WI  = EES_WI(:,:,all(~isnan(EES_WI),[1 2]));%remove nans trial pairs
            
                EMS_WI = WI.EMS.BT_ES(:,:,logicalIdx2);
                EMS_BI = BI.EMS.BT_ES(:,:,:);
                EMS_WI  = EMS_WI(:,:,all(~isnan(EMS_WI),[1 2]));%remove nans trial pairs
                EMS_BI  = EMS_BI(:,:,all(~isnan(EMS_BI),[1 2]));%remove nans trial pairs
            
            
                EES_diff = calc_diff_all_pairs(EES_WI, EES_BI, 10000);
                EMS_diff = calc_diff_all_pairs(EMS_WI, EMS_BI, 10000);
            
                [EES_t_values, EES_p_values] = calc_ttest(EES_diff);
                [EMS_t_values, EMS_p_values] = calc_ttest(EMS_diff);
            
                [EES_surrogate_t_values, EES_surrogate_p_values] = calc_surrogate(EES_WI, EES_BI, 1000);
            
                [surrogate_largest_cluster_sizes, surrogate_largest_cluster_sums] = getLargestClusterSizes(EES_surrogate_p_values<0.05, EES_surrogate_t_values);
                [real_largest_cluster_sizes, real_largest_cluster_sums] = getLargestClusterSizes(EES_p_values<0.05, EES_t_values);
                % change from
            
                final_p_by_channels_image(chan_idx, image_id) = (sum( surrogate_largest_cluster_sums < real_largest_cluster_sums) / length(surrogate_largest_cluster_sums)) * 100;
                % anat_by_channels(chan_idx) = unique(string(BI.EES.chan_ctrl_table.anatomical_label));
            
            
            end
        end
    
        plot_p_values(final_p_by_channels, anat_by_channels)
        sgtitle([sprintf('%s Patient: %d, Enc: %d, Image: %s', ...
            params.btwn_trial_type, patient_id, enc_id, C{1, image_id})]);
        save(sprintf("p_values_by_channel_%s.mat", num2str(patient_id)), "final_p_by_channels", "anat_by_channels", "channels_ids_to_use")
        % save(sprintf("p_values_by_channel_%s_%s.mat", num2str(patient_id), num2str(image_id)), "final_p_by_channels", "anat_by_channels", "channels_ids_to_use")
        save_p_values(patient_id, comp_options, enc_id, image_id, params, final_p_by_channels, anat_by_channels, channel_ids_to_use)
    end
end

% moved to subfunctions
% function plot_p_values(final_p_by_channels, anat_by_channels)
%     figure;
%     % Sort by descending order of final_p_by_channels
%     [sorted_p, sort_idx] = sort(final_p_by_channels, 'descend');
%     sorted_anat = anat_by_channels(sort_idx);
% 
%     % Plot the sorted data
%     bar(sorted_p)
%     xticks(1:length(sorted_anat)); % Ensure all ticks are present
%     xticklabels(sorted_anat)
%     ylabel('p-value')
%     xlabel('Channels')
%     title('Sorted p-values by Channel')
% 
%     % Rotate x-axis labels for better visibility
%     xtickangle(90)
% end

function [new_WI_data, new_BI_data] = shuffle_labels(WI_data, BI_data)
    % first and second dimensions of inputs should be t he same.
    % shuffle slices along the 3rd dimension and place them in new
    % matrices
    combined_data = cat(3, WI_data, BI_data); % Concatenate WI_data and BI_data along the 3rd dimension (TRIAL PAIRS)

    % Shuffle slices along the 3rd dimension
    shuffled_indices = randperm(size(combined_data, 3));
    shuffled_data = combined_data(:, :, shuffled_indices);

    % Split shuffled_data back into new_WI_data and new_BI_data
    nWI = size(WI_data, 3); % Number of slices in WI_data
    new_WI_data = shuffled_data(:, :, 1:nWI);
    new_BI_data = shuffled_data(:, :, nWI+1:end);
end

function [all_t_values, all_p_values] = calc_surrogate(WI_data, BI_data, n_perm)
    sanity_check_sizes(WI_data, BI_data)

    all_t_values = nan(size(WI_data, 1), size(WI_data, 2), n_perm);
    all_p_values = nan(size(WI_data, 1), size(WI_data, 2), n_perm);
    for perm = 1:n_perm
        % shuffle WI, BI labels and calculate ttest nperm times
    
        [new_WI_data, new_BI_data] = shuffle_labels(WI_data, BI_data);
        % run t test
        diff = calc_diff_all_pairs(new_WI_data, new_BI_data, 10000);
        [t_values, p_values] = calc_ttest(diff);
        all_t_values(:,:,perm) = t_values;
        all_p_values(:,:,perm) = p_values;
    end
end

function sanity_check_sizes(WI_data, BI_data)
    % make sure that first and second dimensions have same size
    assert(isequal(size(WI_data, 1), size(BI_data, 1)), 'First dimension of inputs must match.');
    assert(isequal(size(WI_data, 2), size(BI_data, 2)), 'Second dimension of inputs must match.');
end

function [largest_cluster_sizes, largest_cluster_sums] = getLargestClusterSizes(matrix3D, Tmatrix3D)
% getLargestClusterSizes - Computes the largest cluster size and sum of T values for each slice along the 3rd dimension.
% 
% Syntax:
%   [largest_cluster_sizes, largest_cluster_sums] = getLargestClusterSizes(matrix3D, Tmatrix3D)
% 
% Inputs:
%   matrix3D - A 3D matrix of 1's and 0's of size (:,:,N).
%   Tmatrix3D - A 3D matrix of corresponding T values of the same size as matrix3D.
% 
% Outputs:
%   largest_cluster_sizes - A vector of length N, where each element 
%                           corresponds to the size of the largest cluster 
%                           in the respective slice of matrix3D.
%   largest_cluster_sums - A vector of length N, where each element 
%                          corresponds to the sum of the T values within the
%                          largest cluster in the respective slice of matrix3D.

    % Get the number of slices along the 3rd dimension
    numSlices = size(matrix3D, 3);

    % Initialize the output vectors
    largest_cluster_sizes = zeros(numSlices, 1);
    largest_cluster_sums = zeros(numSlices, 1);

    for k = 1:numSlices
        % Extract the current 2D slice
        currentSlice = matrix3D(:, :, k);
        currentT = Tmatrix3D(:, :, k);

        % Identify connected components (clusters) in the slice
        connectedComponents = bwconncomp(currentSlice);

        % Compute the sizes of all clusters
        clusterSizes = cellfun(@numel, connectedComponents.PixelIdxList);

        % Find the largest cluster size (or 0 if no clusters exist)
        if ~isempty(clusterSizes)
            [largestSize, largestIdx] = max(clusterSizes);
            largest_cluster_sizes(k) = largestSize;

            % Sum the T values for the largest cluster
            largestClusterPixels = connectedComponents.PixelIdxList{largestIdx};
            largest_cluster_sums(k) = sum(currentT(largestClusterPixels));
        else
            largest_cluster_sizes(k) = 0;
            largest_cluster_sums(k) = 0;
        end
    end
end

function channel_IDs = load_channels_from_file(patient_id, comp_options, enc_id, image_id, params)
    file_to_load = sprintf('%s\\%s\\%s\\%s\\enc%s_image%s\\all_channels_data.mat', ...
        params.output_folder, num2str(patient_id), comp_options{1}, params.btwn_trial_type, num2str(enc_id), num2str(image_id));
    if exist(file_to_load, 'file')
        channel_IDs = load(file_to_load, 'unique_channel_IDs').unique_channel_IDs;
    else
        sprintf("missing all_channels_data.mat for %s", num2str(patient_id))
        channel_IDs = [];
    end
end

function channel_ids_to_keep = filter_channels_ids_to_use(channel_ids_to_use, original_labels)
    % Initialize the output array
    channel_ids_to_keep = [];
    
    % List of terms to exclude
    excluded_terms = ["pallidum", "putamen", "caudate", "ventricle", "white matter", "insula"];
    
    % Loop through the channel IDs
    for chan_id = channel_ids_to_use
        % Convert the label to lowercase for case-insensitive comparison
        label = lower(string(original_labels(chan_id)));
        
        % Check if the label does not contain any of the excluded terms
        if all(~contains(label, excluded_terms))
            % Append the channel ID to the output array
            channel_ids_to_keep = [channel_ids_to_keep, chan_id];
        end
    end
end
