%% desription of the test:
% have within_image trial enc-delay similarities & between_image. 
% Pick 10000 WI vs BI pairs and compute the difference,
%v4 parallelize over channel*image combinations
%v3 keep just the across-patient & image processing 
%v2 adding image_id dimension to t_test analysis

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
% params.patient_IDs
% if patient_id == 201901
%     channel_ids_to_use = [10 11 12 22 61 7 8 9 22 23 24 25 26 48 49 50 63 64 65 66];
% elseif patient_id == 201908
%     channel_ids_to_use = [28, 30, 37, 38, 74, 76, 82, 83, 84, 85];
% elseif patient_id == 201910
%     channel_ids_to_use = [7 14  25 26 27  41 42  52 53 54  59 60];
% end

addpath("C:\Users\drfrbc\OneDrive - University of Missouri\data\RSA_analysis\Code\Power Spectrum Similarity\subfunctions")
% Comparison options and titles

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

function save_p_values(patient_id, enc_id, p_values_by_chan, anat_by_chan, chan_id_by_chan)
    % save(sprintf('%s\\%s\\%s\\%s\\enc%s_image%s\\t_test_p_values.mat',  params.output_folder, num2str(patient_id), comp_options{2}, params.btwn_trial_type, num2str(enc_id), num2str(image_id)), ...
    %     'p_values_by_chan', 'anat_by_chan', 'chan_id_by_chan')
    save(sprintf("results/t_testing/%s_%s.mat", patient_id, enc_id), "p_values_by_chan", "anat_by_chan", "chan_id_by_chan")
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
comp_options = {'corr BT BI', 'corr BT WI'};
comp_titles = {'Between Trials - Between Images', 'Between Trials - Within Images'};
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

%% do significance testing using 'surrogate'
comp_options = {'corr BT BI', 'corr BT WI'};
enc_id = 1;
% for btwn_trial_type = {'EMS'}%,'EES'}
% params.btwn_trial_type = btwn_trial_type{1};
for patient_id = params.patient_IDs
    fprintf("computing p values for p%s\n", num2str(patient_id))

    % load image_labels, anat_labels
    patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
    % image_labels = load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C');
    anat_labels = load(fullfile(patient_preprocessed_data_path, "D_OWM_t_bipolar.mat"), 'labelsanatbkedit');
    clear patient_preprocessed_data_path

    % load channels_ids_to_use, corresponding anat_labels
    channel_ids_to_use = extractIntegersFromFilenames(patient_id, comp_options, enc_id, 1, params); % image_id=1 is used as a default here. Channels should change with respect to iamge_id.
    original_anat_labels = anat_labels.labelsanatbkedit.anatmacro1(channel_ids_to_use);
    channel_ids_to_use = filter_channels_ids_to_use(channel_ids_to_use, original_anat_labels);
    anat_by_channels =  anat_labels.labelsanatbkedit.anatmacro1(channel_ids_to_use);

    num_images = 9;
    num_channels = length(channel_ids_to_use);
    
    % Preallocate the structure to store results
    final_p_by_channels_image = struct();
    final_p_by_channels_image.EMS = nan(num_channels, num_images);
    
    % Create a combined index for parfor
    num_combinations = num_images * num_channels;
    
    % for comb_idx = 1:num_combinations
        % Decode the combination index into image and channel indices
        % [image_id, chan_idx] = ind2sub([num_images, num_channels], comb_idx);
    
        % Get the channel ID
        % chan_id = channel_ids_to_use(chan_idx);

        chan_id = 20;
        image_id = 2;
    
        % Load BT-BI data for this channel
        WI = struct();
        BI = struct();
        [BI.EMS, WI.EMS] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params, 'EMS');
    
        % Create logical indices to omit same-trial pairs for WI
        [logicalIdx1, logicalIdx2] = createLogicalIndex(size(WI.EMS.BT_ES, 3));
    
        % Process EMS data
        EMS_WI = WI.EMS.BT_ES(:, :, logicalIdx2);
        EMS_BI = BI.EMS.BT_ES(:, :, :);
        EMS_WI = EMS_WI(:, :, all(~isnan(EMS_WI), [1 2])); % Remove NaNs from WI
        EMS_BI = EMS_BI(:, :, all(~isnan(EMS_BI), [1 2])); % Remove NaNs from BI
    
        fprintf("Computing p-values for patient %d, image %d, channel %d\n", patient_id, image_id, chan_id);

        %%
        figure;
        subplot(1,2,1)
        imagesc(mean(EMS_WI,3))
        clim([-0.15, 0.3])
        subplot(1,2,2)
        imagesc(mean(EMS_BI,3))
        clim([-0.15, 0.3])
        sgtitle(sprintf("image%s chan%s %s", num2str(image_id), num2str(chan_id), anat_labels.labelsanatbkedit.anatmacro1(chan_id)))
        %%
    
        % Compute p-values and store in the preallocated structure
        EMS_results(comb_idx) = do_t_testing(EMS_WI, EMS_BI);
    % end
    % Translate EMS_results back to a 2D array after parfor
    final_p_by_channels_image = struct();
    final_p_by_channels_image.EMS = reshape(EMS_results, num_channels, num_images);

    save_p_values(patient_id, enc_id, final_p_by_channels_images, anat_by_channels, channel_ids_to_use)
    print("saved p values for %s", num2str(patient_id))
end
% end

function final_p = do_t_testing(WI, BI)
    diff = calc_diff_all_pairs(WI, BI, 10000);
    [t_values, p_values] = calc_ttest(diff);
    [surrogate_t_values, surrogate_p_values] = calc_surrogate(WI, BI, 1000);
    % figure;
    % subplot(2,2,1)
    % imagesc(p_values)
    % subplot(2,2,2)
    % imagesc(p_values<0.05)
    % subplot(2,2,3)
    % imagesc(mean(surrogate_p_values,3))    
    % subplot(2,2,4)
    % imagesc(mean(surrogate_p_values,3)<0.05)
    [surrogate_largest_cluster_sizes, surrogate_largest_cluster_sums] = getLargestClusterSizes(surrogate_p_values<0.05, surrogate_t_values);
    [real_largest_cluster_sizes, real_largest_cluster_sums] = getLargestClusterSizes(p_values<0.05, t_values);
    final_p = (sum( surrogate_largest_cluster_sums < real_largest_cluster_sums) / length(surrogate_largest_cluster_sums)) * 100;
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
            largestClusterPixels = connectedComponents.PixelIdxList{largestIdx}; % TODO: plot the largest cluster
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
    for chan_idx = 1:length(channel_ids_to_use)
        chan_id = channel_ids_to_use(chan_idx);
        % Convert the label to lowercase for case-insensitive comparison
        label = lower(string(original_labels(chan_idx)));
        
        % Check if the label does not contain any of the excluded terms
        if all(~contains(label, excluded_terms))
            % Append the channel ID to the output array
            channel_ids_to_keep = [channel_ids_to_keep, chan_id];
        end
    end
end