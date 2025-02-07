addpath("../subfunctions")
params = get_parameters();
type = 'EMS';
average_diff = true; % false will do individual trial pairs
comp_options = {'corr BT BI', 'corr BT WI'}; % used to load BI WI.
enc_id = 1;
image_id = 1;
patient_ids = [201901, 201908, 201915, 201906, 201903];
% patient_ids = [201915];



combined_results_table = table();
% for row_id = 1%:10
%     patient_ids = [selected_area_table.patient_id(row_id)];
%     channel_ids_to_use = [selected_area_table.chan_id(row_id)];
%     image_id = selected_area_table.image_id(row_id);
for i = 1:length(patient_ids)
    patient_id = patient_ids(i);
    for image_id = 1:9
        
    patient_preprocessed_data_path = fullfile(params.preprocessed_data_location, sprintf('/CS%s/', num2str(patient_id)));
    image_labels = load(fullfile(patient_preprocessed_data_path, "OWM_trialinfo.mat"), 'C');
    anat_labels = load(fullfile(patient_preprocessed_data_path, "D_OWM_t_bipolar.mat"), 'labelsanatbkedit');
    
    channel_ids_to_use = extractIntegersFromFilenames(patient_id, comp_options, enc_id, image_id, params); % image_id=1 is used as a default here. Channels should change with respect to iamge_id.
    original_anat_labels = anat_labels.labelsanatbkedit.anatmacro1(channel_ids_to_use);
    channel_ids_to_use = filter_channels_ids_to_use(channel_ids_to_use, original_anat_labels);
    % channel_ids_to_use = [65, 28, 13, 57, 14]
    anat_by_channels =  anat_labels.labelsanatbkedit.anatmacro1(channel_ids_to_use);
    p_by_channels = nan(length(channel_ids_to_use), 1);
    
    for chan_idx = 1:length(channel_ids_to_use)
        chan_id = channel_ids_to_use(chan_idx);

% chan_id = 20;
% anat = string(anat_labels.labelsanatbkedit.anatmacro1(chan_id));
    
        [WI, BI] = get_and_pre_process_WI_BI(patient_id, comp_options, enc_id, image_id, chan_id, params, type);
        WI=mean(WI, [1 2]);
        BI=mean(BI, [1 2]);
        if average_diff
            diff = calc_diff_avg(WI, BI, 1000);
        else
            diff = calc_diff_all_pairs(WI, BI, 10000);
        end
        [real_t_values, real_p_values] = calc_ttest(diff);
        [surrogate_t_values, surrogate_p_values] = calc_surrogate(WI, BI, 1000, average_diff);
    
        % final_p = cluster_analysis_and_visualization_abs(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values);
        final_p = mean(real_t_values > surrogate_t_values);
        p_by_channels(chan_idx) = final_p;
    end

    chan_ids = channel_ids_to_use;
    p_values = p_by_channels;
    anat_labels_by_channels = anat_by_channels;
    
    % Create the table
    result_table = table(repmat(patient_id, length(channel_ids_to_use), 1), ...
                         channel_ids_to_use', p_by_channels, anat_by_channels, ...
                         repmat(image_id, length(channel_ids_to_use), 1), ...
                         'VariableNames', {'patient_id', 'chan_id', 'p', 'anat', 'image_id'});
    
    result_table = sortrows(result_table, 'p','descend');

    combined_results_table = [combined_results_table; result_table];
    end

end
%%
 combined_results_table.patient_id = string(combined_results_table.patient_id);
combined_results_table.anat = string(combined_results_table.anat);
 combined_results_table = sortrows(combined_results_table,'p','descend');

summary_table = groupsummary(combined_results_table, {'patient_id', 'chan_id', 'anat'}, ...
{'mean', 'std'}, 'p');

summary_table = sortrows(summary_table,'mean_p','descend');
save("summary_p_tables_by_t.mat", "summary_table", "combined_results_table")

% % Get unique anat labels for each (patient_id, chan_id)
% anat_labels = rows2vars(combined_results_table(:, {'patient_id', 'chan_id', 'anat'}));
% anat_labels_unique = unique(anat_labels, 'rows', 'stable'); % Make sure the order is preserved
% 
% % Check if there is only one unique anat label per (patient_id, chan_id)
% disp(anat_labels_unique);
% 
% % Merge anat labels with the summary_table
% summary_table_with_anat = join(summary_table, anat_labels_unique, ...
%                                'Keys', {'patient_id', 'chan_id'}, ...
%                                'RightVariables', {'anat'});
% 
% % Display the updated table
% disp(summary_table_with_anat);
%%
% final_p = mean(real_p_values > surrogate_p_values);

%%

% testing the line below
% cluster_analysis_and_visualization(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values)
% cluster_analysis_and_visualization_abs(real_t_values, real_p_values, surrogate_t_values, surrogate_p_values)

surrogate_clusters = struct();
surrogate_clusters.positive_t = struct();
surrogate_clusters.negative_t = struct();
[surrogate_clusters.positive_t.largest_cluster_sizes, surrogate_clusters.positive_t.largest_cluster_sums, surrogate_clusters.positive_t.largest_cluster_locations] = getLargestClusterSizes(surrogate_p_values<0.05&surrogate_t_values>0, surrogate_t_values);
[surrogate_clusters.negative_t.largest_cluster_sizes, surrogate_clusters.negative_t.largest_cluster_sums, surrogate_clusters.negative_t.largest_cluster_locations] = getLargestClusterSizes(surrogate_p_values<0.05&surrogate_t_values<0, surrogate_t_values);

real_clusters = struct();
real_clusters.positive_t = struct();
real_clusters.negative_t = struct();
[real_clusters.positive_t.largest_cluster_sizes, real_clusters.positive_t.largest_cluster_sums, real_clusters.positive_t.largest_cluster_locations] = getLargestClusterSizes(real_p_values<0.05&real_t_values>0, real_t_values);
[real_clusters.negative_t.largest_cluster_sizes, real_clusters.negative_t.largest_cluster_sums, real_clusters.negative_t.largest_cluster_locations] = getLargestClusterSizes(real_p_values<0.05&real_t_values<0, real_t_values);

% [real_cluster_sizes, real_cluster_sums, real_cluster_locations] = getAllClusterData(real_p_values<0.05, real_t_values);

% final_p = struct();
% final_p.positive_t = (sum( surrogate_clusters.positive_t.largest_cluster_sums > real_clusters.positive_t.largest_cluster_sums) / length(surrogate_clusters.positive_t.largest_cluster_sums)) * 100;
% final_p.negative_t = (sum( surrogate_clusters.negative_t.largest_cluster_sums < real_clusters.negative_t.largest_cluster_sums) / length(surrogate_clusters.negative_t.largest_cluster_sums)) * 100;
% % final_p = (sum( surrogate_largest_cluster_sums < real_largest_cluster_sums) / length(surrogate_largest_cluster_sums)) * 100;

% save_name = sprintf("p_data_p%s_image%s_chan%s",num2str(patient_id),num2str(image_id),num2str(chan_id));
% test=0;
% while exist(sprintf('%s_test%s.mat',save_name, num2str(test)), 'file')
%     test = test + 1;
% end
% save_name = sprintf('%s_test%s.mat',save_name, num2str(test));
save_name = sprintf("results/selected_p_data_p%s_image%s_chan%s.mat",num2str(patient_id),num2str(image_id),num2str(chan_id));
save(save_name, "final_p", "real_t_values", "real_p_values", "real_clusters","surrogate_t_values", "surrogate_p_values", "surrogate_clusters")
%% testing; get all clusters for the real data
% Modify the cluster processing step to include both all and largest clusters
all_clusters = struct();
largest_clusters = struct();

% Compute clusters for positive t-values
[largest_clusters.positive_t.largest_cluster_sizes, ...
 largest_clusters.positive_t.largest_cluster_sums, ...
 largest_clusters.positive_t.largest_cluster_locations] = ...
    getLargestClusterSizes(real_p_values < 0.05 & real_t_values > 0, real_t_values);

[all_clusters.positive_t.all_cluster_sizes, ...
 all_clusters.positive_t.all_cluster_sums, ...
 all_clusters.positive_t.all_cluster_locations] = ...
    getAllClusterData(real_p_values < 0.05 & real_t_values > 0, real_t_values);

% Compute clusters for negative t-values
[largest_clusters.negative_t.largest_cluster_sizes, ...
 largest_clusters.negative_t.largest_cluster_sums, ...
 largest_clusters.negative_t.largest_cluster_locations] = ...
    getLargestClusterSizes(real_p_values < 0.05 & real_t_values < 0, real_t_values);

[all_clusters.negative_t.all_cluster_sizes, ...
 all_clusters.negative_t.all_cluster_sums, ...
 all_clusters.negative_t.all_cluster_locations] = ...
    getAllClusterData(real_p_values < 0.05 & real_t_values < 0, real_t_values);

% Threshold for largest clusters
positive_threshold = prctile(surrogate_clusters.positive_t.largest_cluster_sums, 95);
negative_threshold = prctile(surrogate_clusters.negative_t.largest_cluster_sums, 95);

% Visualize positive t clusters
visualize_significant_clusters_WI_vs_BI_2(real_t_values, ...
    all_clusters.positive_t, ...
    largest_clusters.positive_t, ...
    positive_threshold, alpha);

% Visualize negative t clusters
visualize_significant_clusters_WI_vs_BI_2(real_t_values, ...
    all_clusters.negative_t, ...
    largest_clusters.negative_t, ...
    negative_threshold, alpha);
%% plotting
plot_similarity_means(WI,BI)
sgtitle(sprintf("image%s chan%s %s", num2str(image_id), num2str(chan_id), anat))

plot_t_with_p_highlights(real_t_values, real_p_values<0.05|real_t_values>0)

% plot_p_values(real_p_values, surrogate_p_values) % plot
plot_t_p_values(real_p_values, surrogate_p_values, real_t_values, surrogate_t_values)
sgtitle(sprintf("image%s chan%s %s", num2str(image_id), num2str(chan_id), anat))

% plot size pdf
% positive t
plot_cluster_size_pdf(real_clusters.positive_t.largest_cluster_sums, surrogate_clusters.positive_t.largest_cluster_sums, 1000);
sgtitle(sprintf("positive t image%s chan%s %s", num2str(image_id), num2str(chan_id), anat))
% negative t
plot_cluster_size_pdf(real_clusters.negative_t.largest_cluster_sums, surrogate_clusters.negative_t.largest_cluster_sums, 1000);
sgtitle(sprintf("negative t image%s chan%s %s", num2str(image_id), num2str(chan_id), anat))

% overlay significant clusters onto
alpha = 0.05;   % Significance level
% threshold = prctile(real_largest_cluster_sizes, 95);
% visualize_significant_clusters(real_p_values, real_largest_cluster_sizes, threshold, alpha);

% threshold = prctile(surrogate_largest_cluster_sums, 95);
% visualize_significant_clusters_WI_vs_BI(real_t_values, real_largest_cluster_sums, threshold, alpha);

%positive t
threshold = prctile(surrogate_clusters.positive_t.largest_cluster_sums, 95);
visualize_significant_clusters_WI_vs_BI(real_t_values, real_clusters.positive_t.largest_cluster_sums, threshold, alpha);
sgtitle(sprintf("image%s chan%s %s", num2str(image_id), num2str(chan_id), anat))
% negative t
threshold = prctile(surrogate_clusters.negative_t.largest_cluster_sums, 95);
visualize_significant_clusters_WI_vs_BI(real_t_values, real_clusters.negative_t.largest_cluster_sums, threshold, alpha);
sgtitle(sprintf("image%s chan%s %s", num2str(image_id), num2str(chan_id), anat))

%% loading
% function [EMS_WI,EMS_BI] = get_and_pre_process_WI_BI(patient_id, comp_options, enc_id, image_id, chan_id, params, type)
%         [BI, WI] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params, type); % doesn't matter that it looks like type='EMS' is returned.
% 
%         % Create logical indices to omit same-trial pairs for WI
%         [logicalIdx1, logicalIdx2] = createLogicalIndex(size(WI.BT_ES, 3));
% 
%         % Process EMS data
%         EMS_WI = WI.BT_ES(:, :, logicalIdx2);
%         EMS_BI = BI.BT_ES(:, :, :);
%         EMS_WI = EMS_WI(:, :, all(~isnan(EMS_WI), [1 2])); % Remove NaNs from WI
%         EMS_BI = EMS_BI(:, :, all(~isnan(EMS_BI), [1 2])); % Remove NaNs from BI
% end

%% computations
% function diff = calc_diff_all_pairs(WI, BI, num_pairs)
%     % Inputs:
%     % WI - The within-item matrix (size: time x time x num_WI)
%     % BI - The between-item matrix (size: time x time x num_BI)
%     % num_pairs - Maximum number of pairs to compute (randomly selected if less than total_pairs)
% 
%     % Dimensions of the input matrices
%     num_WI = size(WI, 3);
%     num_BI = size(BI, 3);
% 
%     % Determine the total number of pairs
%     total_pairs = num_WI * num_BI;
% 
%     % Validate num_pairs
%     if nargin < 3 || isempty(num_pairs)
%         num_pairs = total_pairs; % Default to all pairs if num_pairs not specified
%     end
%     num_pairs = min(num_pairs, total_pairs);
% 
%     % Generate random pair indices if num_pairs is less than total_pairs
%     if num_pairs < total_pairs
%         % rng(0); % For reproducibility (optional, remove if not needed)
%         random_indices = randperm(total_pairs, num_pairs);
%     else
%         random_indices = 1:total_pairs; % Use all combinations if num_pairs equals total_pairs
%     end
% 
%     % Preallocate the output difference matrix
%     diff = zeros(size(WI, 1), size(WI, 2), num_pairs);
% 
%     % Calculate differences for the selected combinations
%     for k = 1:num_pairs
%         [i, j] = ind2sub([num_WI, num_BI], random_indices(k)); % Map linear index to subscripts
%         diff(:, :, k) = WI(:, :, i) - BI(:, :, j);
%     end
% end

%% computations
% function diff = calc_diff_avg(WI, BI, num_pairs)
%     % Inputs:
%     % WI - The within-item matrix (size: time x time x num_WI)
%     % BI - The between-item matrix (size: time x time x num_BI)
%     % num_pairs - number of pairs to compute (randomly selected if less than total_pairs)
% 
%     % Dimensions of the input matrices
%     num_WI = size(WI, 3);
%     num_BI = size(BI, 3);
% 
%     % Calculate the minimum size for the 3rd dimension
%     min_size = min(num_WI, num_BI);
% 
%     % Generate random sets of indices to select from the 3rd dimension
%     if num_WI < num_BI
%         % Generate `num_pairs` sets of random indices for WI (1:num_WI) and BI (min_size random from 1:num_BI)
%         random_indices_WI = arrayfun(@(x) 1:num_WI, 1:num_pairs, 'UniformOutput', false);
%         random_indices_BI = arrayfun(@(x) randperm(num_BI, min_size), 1:num_pairs, 'UniformOutput', false);
%     else
%         % Generate `num_pairs` sets of random indices for WI (min_size random from 1:num_WI) and BI (1:num_BI)
%         random_indices_WI = arrayfun(@(x) randperm(num_WI, min_size), 1:num_pairs, 'UniformOutput', false);
%         random_indices_BI = arrayfun(@(x) 1:num_BI, 1:num_pairs, 'UniformOutput', false);
%     end
% 
%     % this method chooses random, repeating
%     % random_indices_WI = randi(num_WI, [num_pairs, min_size]);  % num_pairs sets of random indices for WI
%     % random_indices_BI = randi(num_BI, [num_pairs, min_size]);  % num_pairs sets of random indices for BI
% 
%     % Preallocate the output difference matrix
%     diff = zeros(size(WI, 1), size(WI, 2), num_pairs);
% 
%     % Calculate differences for the selected combinations
%     for k = 1:num_pairs
%         % Select the random indices for WI and BI
%         % indices_WI = random_indices_WI(k, :); % for randi, repeating
%         % indices_BI = random_indices_BI(k, :);
% 
%         indices_WI = random_indices_WI{k}; % for randperm
%         indices_BI = random_indices_BI{k};
% 
%         % Compute mean across the selected indices for WI and BI
%         WI_mean = mean(WI(:, :, indices_WI), 3);
%         BI_mean = mean(BI(:, :, indices_BI), 3);
% 
%         % Compute the difference
%         diff(:, :, k) = WI_mean - BI_mean;
%     end
% end


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

% 
% function [all_t_values, all_p_values, all_diffs] = calc_surrogate(WI_data, BI_data, n_perm, average_diff)
%     sanity_check_sizes(WI_data, BI_data)
% 
%     if average_diff
%         all_diffs = nan(nperm, 1000);
%     else
%         all_diffs = nan(nperm, 10000);
%     end
% 
%     all_t_values = nan(size(WI_data, 1), size(WI_data, 2), n_perm);
%     all_p_values = nan(size(WI_data, 1), size(WI_data, 2), n_perm);
%     for perm = 1:n_perm
%         % shuffle WI, BI labels and calculate ttest nperm times
%         [new_WI_data, new_BI_data] = shuffle_labels(WI_data, BI_data);
% 
%         % calculate differences
%         if average_diff
%             diff = calc_diff_avg(new_WI_data, new_BI_data, 1000);
%         else
%             diff = calc_diff_all_pairs(new_WI_data, new_BI_data, 10000);
%         end
%         all_diffs(perm, :) = diff;
% 
%         % run t test between diff and 0
%         [t_values, p_values] = calc_ttest(diff);
%         all_t_values(:,:,perm) = t_values;
%         all_p_values(:,:,perm) = p_values;
%     end
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

function [largest_cluster_sizes, largest_cluster_sums, largest_cluster_locations] = getLargestClusterSizes(matrix3D, Tmatrix3D)
% getLargestClusterSizes - Computes the largest cluster size, sum of T values, 
% and locations of the largest cluster for each slice along the 3rd dimension.
% 
% Syntax:
%   [largest_cluster_sizes, largest_cluster_sums, largest_cluster_locations] = getLargestClusterSizes(matrix3D, Tmatrix3D)
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
%   largest_cluster_locations - A 3D logical matrix where each slice (:,:,k) 
%                               indicates the pixels belonging to the largest cluster 
%                               in the respective slice of matrix3D.

    % Get the number of slices along the 3rd dimension
    numSlices = size(matrix3D, 3);

    % Initialize the output variables
    largest_cluster_sizes = zeros(numSlices, 1);
    largest_cluster_sums = zeros(numSlices, 1);
    largest_cluster_locations = false(size(matrix3D));

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

            % Get the indices of the largest cluster
            largestClusterPixels = connectedComponents.PixelIdxList{largestIdx};

            % Create a logical matrix for the largest cluster
            largestClusterMask = false(size(currentSlice));
            largestClusterMask(largestClusterPixels) = true;
            largest_cluster_locations(:, :, k) = largestClusterMask;

            % Sum the T values for the largest cluster
            largest_cluster_sums(k) = sum(currentT(largestClusterPixels));
        else
            largest_cluster_sizes(k) = 0;
            largest_cluster_sums(k) = 0;
        end
    end
end

function [all_cluster_sizes, all_cluster_sums, all_cluster_locations] = getAllClusterData(matrix3D, Tmatrix3D)
% getAllClusterData - Computes all cluster sizes, sums of T values, 
% and locations of all clusters for each slice along the 3rd dimension.
% 
% Syntax:
%   [all_cluster_sizes, all_cluster_sums, all_cluster_locations] = getAllClusterData(matrix3D, Tmatrix3D)
% 
% Inputs:
%   matrix3D - A 3D matrix of 1's and 0's of size (:,:,N).
%   Tmatrix3D - A 3D matrix of corresponding T values of the same size as matrix3D.
% 
% Outputs:
%   all_cluster_sizes - A cell array of length N, where each cell contains a vector 
%                       of sizes of all clusters in the respective slice of matrix3D.
%   all_cluster_sums - A cell array of length N, where each cell contains a vector 
%                      of the sums of the T values for all clusters in the respective slice of matrix3D.
%   all_cluster_locations - A cell array of length N, where each cell contains a cell array 
%                           of indices (pixel locations) for all clusters in the respective slice of matrix3D.

    % Get the number of slices along the 3rd dimension
    numSlices = size(matrix3D, 3);

    % Initialize the output variables
    all_cluster_sizes = cell(numSlices, 1);
    all_cluster_sums = cell(numSlices, 1);
    all_cluster_locations = cell(numSlices, 1);

    for k = 1:numSlices
        % Extract the current 2D slice
        currentSlice = matrix3D(:, :, k);
        currentT = Tmatrix3D(:, :, k);

        % Identify connected components (clusters) in the slice
        connectedComponents = bwconncomp(currentSlice);

        % Compute the sizes of all clusters
        clusterSizes = cellfun(@numel, connectedComponents.PixelIdxList);

        % Initialize cluster sums and locations
        clusterSums = zeros(length(clusterSizes), 1);
        clusterLocations = cell(length(clusterSizes), 1);

        % Loop through each cluster to compute sums and store locations
        for i = 1:length(clusterSizes)
            clusterPixels = connectedComponents.PixelIdxList{i};
            clusterSums(i) = sum(currentT(clusterPixels));
            clusterLocations{i} = clusterPixels;
        end

        % Store results for the current slice
        all_cluster_sizes{k} = clusterSizes;
        all_cluster_sums{k} = clusterSums;
        all_cluster_locations{k} = clusterLocations;
    end
end


%% plotting
% function plot_similarity_means(WI,BI)
%         figure;
%         subplot(1,2,1)
%         imagesc(mean(WI,3))
%         clim([-0.15, 0.3])
%         title(sprintf('WI mean %s', num2str(size(WI,3))))
% 
% 
%         subplot(1,2,2)
%         imagesc(mean(BI,3))
%         clim([-0.15, 0.3])
%         title(sprintf('BI mean %s', num2str(size(BI,3))))
% end

function plot_p_values(real_p_values, surrogate_p_values)
    figure;
    subplot(2,2,1)
    imagesc(real_p_values)
    title('real p values')
    clim([0 1])

    subplot(2,2,2)
    imagesc(real_p_values<0.05)
    title('real p values < 0.05')
    clim([0 1])

    subplot(2,2,3)
    imagesc(mean(surrogate_p_values,3))
    title(sprintf('surrogate(shuffled) p values mean (n=%s)', num2str(size(surrogate_p_values,3))))
    clim([0 1])

    subplot(2,2,4)
    imagesc(mean(surrogate_p_values<0.05,3))
    title(sprintf('surrogate(shuffled) p values < 0.05 mean (n=%s)', num2str(size(surrogate_p_values,3))))
    clim([0 1])
end

function plot_t_p_values(real_p_values, surrogate_p_values, real_t_values,  surrogate_t_values)
    figure;
    subplot(3,2,1)
    imagesc(real_t_values)
    title('real t values')
    clim([-20 20])

    subplot(3,2,2)
    imagesc(real_p_values)
    title('real p values')
    clim([0 1])

    subplot(3,2,3)
    imagesc(real_p_values<0.05)
    title('real p values < 0.05')
    clim([0 1])

    subplot(3,2,4)
    imagesc(mean(surrogate_t_values, 3))
    title(sprintf('surrogate(shuffled) t values mean (n=%s)', num2str(size(surrogate_t_values,3))))
    clim([-20 20])

    subplot(3,2,5)
    imagesc(mean(surrogate_p_values,3))
    title(sprintf('surrogate(shuffled) p values mean (n=%s)', num2str(size(surrogate_p_values,3))))
    clim([0 1])

    subplot(3,2,6)
    imagesc(mean(surrogate_p_values<0.05,3))
    title(sprintf('surrogate(shuffled) p values < 0.05 mean (n=%s)', num2str(size(surrogate_p_values,3))))
    clim([0 1])
end

function visualize_significant_clusters_WI_vs_BI(real_t_values, real_cluster_sizes, threshold, alpha)
    figure;
    imagesc(real_t_values);
    hold on;
    colormap('gray');
    h = colorbar;
    ylabel(h, 'P-Value');
    clim([0 alpha]);

    % Find connected components
    real_mask = real_t_values < alpha;
    connectedComponents = bwconncomp(real_mask);
    significant_clusters = real_cluster_sizes > threshold;
    plotted_non_significant = false;
    plotted_significant = false;
    
    % Plot clusters
    for i = 1:length(real_cluster_sizes)
        cluster_indices = connectedComponents.PixelIdxList{i};
        [row, col] = ind2sub(size(real_t_values), cluster_indices);
        
        if significant_clusters(i)
            plot(col, row, 'r.', 'MarkerSize', 10);
            plotted_significant = true;
        else
            plot(col, row, 'b.', 'MarkerSize', 5);
            plotted_non_significant = true;
        end
    end
    
    title('Clusters with Significant T Sum');
    xlabel('Maintenance Indices');
    ylabel('Encoding Indices');
    % Adjust legend based on plotted clusters
    if plotted_significant && plotted_non_significant
        legend('Significant Cluster', 'Non-Significant Cluster', 'Location', 'Best');
    elseif plotted_significant
        legend('Significant Cluster', 'Location', 'Best');
    elseif plotted_non_significant
        legend('Non-Significant Cluster', 'Location', 'Best');
    end
    grid on;
    hold off;
end

function visualize_significant_clusters_WI_vs_BI_2(real_t_values, all_clusters, largest_clusters, threshold, alpha)
    % Visualize all clusters and highlight the largest clusters
    figure;
    imagesc(real_t_values);
    hold on;
    colormap('jet');
    colorbar;
    clim([-20 20]);
    title('Significant Clusters Visualization');
    xlabel('Maintenance Indices');
    ylabel('Encoding Indices');
    
    % Plot all clusters (blue)
    for i = 1:length(all_clusters.all_cluster_sizes)
        cluster_pixels = all_clusters.all_cluster_locations{i};
        [row, col] = ind2sub(size(real_t_values), cluster_pixels);
        plot(col, row, 'b.', 'MarkerSize', 5); % Blue for all clusters
    end
    
    % Highlight largest clusters (red)
    for i = 1:length(largest_clusters.largest_cluster_sizes)
        if largest_clusters.largest_cluster_sums(i) > threshold
            cluster_pixels = largest_clusters.largest_cluster_locations(:, :, i);
            [row, col] = find(cluster_pixels);
            plot(col, row, 'r.', 'MarkerSize', 10); % Red for largest clusters
        end
    end
    
    legend('All Clusters', 'Largest Clusters', 'Location', 'Best');
    hold off;
end

%% helper
% function [logicalIdx1, logicalIdx2] = createLogicalIndex(sizeDim) % omits diagonal for WI pairs
%     % createLogicalIndex creates two logical arrays for indexing a matrix,
%     % excluding diagonal elements (e.g., 1,1; 2,2; ...; n,n).
%     %
%     % INPUT:
%     %   sizeDim - Size of the dimension (e.g., 12 for a matrix (:,:,12,12))
%     %
%     % OUTPUT:
%     %   logicalIdx1 - Logical index for the first dimension
%     %   logicalIdx2 - Logical index for the second dimension
% 
%     % Validate input
%     if ~isscalar(sizeDim) || sizeDim <= 0 || sizeDim ~= round(sizeDim)
%         error('Input sizeDim must be a positive integer.');
%     end
% 
%     % Initialize logical arrays
%     logicalIdx1 = true(sizeDim, sizeDim);
%     logicalIdx2 = true(sizeDim, sizeDim);
% 
%     % Create diagonal indices to exclude
%     diagIdx = 1:sizeDim;
% 
%     % Set diagonal elements to false
%     logicalIdx1(sub2ind([sizeDim, sizeDim], diagIdx, diagIdx)) = false;
%     logicalIdx2(sub2ind([sizeDim, sizeDim], diagIdx, diagIdx)) = false;
% end

function sanity_check_sizes(WI_data, BI_data)
    % make sure that first and second dimensions have same size
    assert(isequal(size(WI_data, 1), size(BI_data, 1)), 'First dimension of inputs must match.');
    assert(isequal(size(WI_data, 2), size(BI_data, 2)), 'Second dimension of inputs must match.');
end