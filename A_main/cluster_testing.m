%%
hellbender = false;
nPermutations = 1000;
use_z = false; % false uses p
alpha = 0.005;%0.05;
z_thresh = 1.96;

plot_cluster_sizes_by_obs = 1:3; % choose obs to plot or set false
plot_sig_clusters_sizes_by_obs = 1:3; % choose obs to plot or set false

%%
if hellbender
    addpath 'Raw Data Storage' %#ok<UNRCH>
    addpath 'subfunctions'
    addpath 'clustering subfunctions'
    output_folder = fullfile('/cluster/VAST/bkybg-lab/Data/OWM Utah Data/RSA/PSS/parallel output/allpatients gammamod allregions allitem allenc');
    patient_IDs = [201907 201908, 201903, 201905, 201906, 201901, 201910, 201915];
else
    addpath('../Z_Raw Data Storage') 
    addpath('../subfunctions')
    addpath('../clustering subfunctions')
    output_folder = 'D:\Power Spectrum Similarity\AA_Processed Data\allpatients gammamod allregions allitem enc1 correct';
    patient_IDs = [201907]; %#ok<NBRAK>
end


%% intialize parallel pool
if isempty(gcp('nocreate')) % Open a local parallel pool if none exists
    if hellbender
        parpool('local', 64) %#ok<UNRCH> % 25-64
    else
        % parpool('local', 2);% %#ok<UNRCH> % 1-2 
    end
end

obs = 1;
for idx = 1:length(patient_IDs)
    patient_ID = patient_IDs(idx);
    PS_file = get_PS_file(output_folder, patient_ID);
    test_file = get_ES_file(PS_file);
    label_table = PS_file.label_table;
    label_table.anatomical_label = string(label_table.anatomical_label);
    rows_without_nan = get_rows_without_nan(label_table);
    % valid_rows = find(rows_without_nan); % indices without nan (orig length)
    label_table = label_table(rows_without_nan,:);
    % % % compute correlations
    for target_enc_ids = 1%:3
        all_matrix = squeeze(test_file.all3_ES_matrix(:,:,:,target_enc_ids));
        all_matrix = all_matrix(rows_without_nan,:,:);
        clear rows_without_nan
        for target_image_ids = 1%:9
            ctrl_file = get_ctrl_file(PS_file, patient_ID, target_enc_ids, target_image_ids);
            unique_channel_IDs = ctrl_file.unique_channel_IDs;
            num_channels = length(unique_channel_IDs);

            % test_rows = valid_rows(ctrl_file.test_rows); % valid_rows corresponds to (all_matrix, 1)
            test_table = label_table(ctrl_file.test_rows,:);
            test_matrix = all_matrix(ctrl_file.test_rows,:,:);
            % test_table = label_table(test_rows,:);

            % control_matrices = ctrl_file.all_control_EMS_matrices_by_chan; % [58x1   cell] nchannelsx1
            clusters_by_anat = struct();
            clusters_by_anat.anat = cell(num_channels);
            clusters_by_anat.num_sig_test_clusters = cell(num_channels);
            clusters_by_anat.sig_test_cluster_sizes = cell(num_channels);
            clusters_by_anat.real_test_cluster_sizes = cell(num_channels);

            clusters_save_file = strrep(ctrl_file.Properties.Source, 'clusters.mat');

            for chan_idx = 1:num_channels
                chan_id = unique_channel_IDs(chan_idx);
                channel_file = matfile(strrep(ctrl_file.Properties.Source, 'all_channels_data.mat', sprintf('BT_%d', chan_id)));
                % control_matrices = ctrl_file.all_control_EMS_matrices_by_chan(chan_idx, 1); % hopefully better slicing than this
                % control_matrices = control_matrices{1}; % for some reason the 4D matrix is still in a cell.
                control_matrices = channel_file.control_EMS_matrices(:, :, :, :); % nEnc, nMaint, ntrials, nSamples (if ntrials==1 then it is averaged already)
                % get test % need to subset by channel?
                % need to combine rows_without_nan with ctrl_file.test_rows
    
                chan_test_rows = filter_by_channel_id(test_table, chan_id);
                chan_test_table = test_table(chan_test_rows,:);
                chan_test_matrix = test_matrix(chan_test_rows,:,:);

                % chan_test_rows = test_rows(chan_test_rows);
                % test_matrix = all_matrix(chan_test_rows,:,:);
                %% now run across trials.
                % null distribution is calculated with ntrialx1 +
                % ntrialx100 choices and the significance should be done on
                % each test trial after.

                % run clustering
                % control: 41 x 640 x nTrialsWItem x 100
                % test: nTrialsWItem x 41 x 640

                % average across TrialsWItem
                control_matrices = squeeze(mean(control_matrices, 3));
                chan_test_matrix = squeeze(mean(chan_test_matrix, 1));
                
                % the test data is w.r.t. 'non_selection_windows' i.e. whole trial except 'selecting
                % object' while the control is w.r.t. maintenance_windows
                if size(control_matrices, 2) ~= size(chan_test_matrix, 2) % probably the case mentioned above
                    % [~, ~, ~, ~, maint_win_IDs, ~, ~] = get_window_IDs(); % won't change each function call
                    % test_matrix = test_matrix(:,maint_win_IDs);
                    % clear maint_win_IDs
                    chan_test_matrix = chan_test_matrix(:,251:640);
                end

                [num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, null_cluster_sizes, threshold_from_null, real_p_val_matrix] = get_clusters(control_matrices, chan_test_matrix, nPermutations, alpha, use_z, z_thresh);

                %% optional plotting
                if any(plot_cluster_sizes_by_obs == obs) % check if this obs should be plotted % test & null cluster size pdfs
                    plot_cluster_size_pdf(real_test_cluster_sizes, null_cluster_sizes, nPermutations); % test & null cluster size PDFs
                    xline(threshold_from_null, 'r', 'LineWidth', 2, 'DisplayName', sprintf('p < %d Threshold from Null', alpha)); % Significance testing
                    clear null_cluster_sizes
                else
                    clear null_cluster_sizes
                end

                if any(plot_sig_clusters_sizes_by_obs == obs) % view test sig/~sig clusters
                    visualize_significant_clusters(real_p_val_matrix, sig_test_cluster_sizes, threshold_from_null, alpha); % view test sig/~sig clusters
                end
                %% gather data across obs to plot by anatomy
                % store {above 3} each iteration & info: anat, channel_ID,
                if length(unique(chan_test_table.channel_ID)) > 1 | length(unique(chan_test_table.anatomical_label)) > 1
                    error("should be 1 anat and 1 chan ID for test data")
                else
                    anat = unique(chan_test_table.anatomical_label);
                end

                % will need to update and use {obs} when computing across
                % encodings.
                clusters_by_anat.anat{chan_idx} = anat;
                clusters_by_anat.num_sig_test_clusters{chan_idx} = num_sig_test_clusters;
                clusters_by_anat.sig_test_cluster_sizes{chan_idx} = sig_test_cluster_sizes;
                clusters_by_anat.real_test_cluster_sizes{chan_idx} = real_test_cluster_sizes;
                
                obs = obs + 1;
            end
            save(clusters_save_file, "clusters_by_anat", "-v7.3")



            %% WIP: plot as a function of anatomy
            [clusters_by_uniq_anat] = sort_by_uniq_anat(clusters_by_anat);

            plot_cluster_sizes_by_uniq_anat(clusters_by_uniq_anat)
            
            
        end
    end
end
%% functions

function rows = filter_by_channel_id(table, target_channel_ID)
    % Function to filter rows based on target channel ID
    % Inputs:
    % - table: the input table containing the field 'channel_ID'
    % - target_channel_ID: the channel ID to filter by

    % Create a logical mask to identify rows with the target channel ID
    rows = table.channel_ID == target_channel_ID;
end

%% updated below to run across trials
% function run_analysis(control_filename, test_filename, nPermutations, alpha)
function [num_sig_test_clusters, sig_test_cluster_sizes, real_test_cluster_sizes, null_cluster_sizes, threshold_from_null, real_p_val_matrix] = get_clusters(control_matrices, test_matrix, nPermutations, alpha, use_z, z_thresh)
    sanity_check_matrices(test_matrix, control_matrices)

    % Compute real test's z-scored and p-value matrices
    % [real_z_scored_mat, real_p_val_matrix] = compute_z_p_matrices(test_matrix, control_matrices);
    [real_z_scored_mat, real_p_val_matrix] = compute_z_p_matrices(test_matrix, control_matrices);
    
    % Compute real test cluster sizes
    if use_z
        real_test_cluster_sizes = compute_cluster_sizes_z(real_z_scored_mat, z_thresh);
    else
        real_test_cluster_sizes = compute_cluster_sizes(real_p_val_matrix, alpha);
    end
    
    % Concatenate control and test matrices for nulldistribution permutations
    all_avg_matrices = cat(3, control_matrices, test_matrix);
    
    % Run permutations
    null_cluster_sizes = run_permutations(nPermutations, all_avg_matrices, alpha, use_z, z_thresh);
    
    % Concatenate null cluster sizes
    null_cluster_sizes = concatenate_cluster_sizes(null_cluster_sizes);

    % find significant
    threshold_from_null = prctile(null_cluster_sizes, 100 * (1 - alpha));
    sig_test_cluster_sizes = real_test_cluster_sizes(real_test_cluster_sizes > threshold_from_null);
    num_sig_test_clusters = numel(sig_test_cluster_sizes); % Get the number of significant clusters
    
end

% function [control_matrices, test_matrix] = load_data(control_filename, test_filename)
    % control_data = load(control_filename);
    % test_data = load(test_filename);
    % 
    % % Extract the matrices
    % control_matrices = control_data.EMS_avg_all_samples;
    % test_matrix = test_data.EMS_avg;

function [z_scored_mat, p_val_matrix] = compute_z_p_matrices(test_matrix, control_matrices)
    [z_scored_mat, p_val_matrix] = get_z_p_mat(test_matrix, control_matrices);
end

function cluster_sizes = compute_cluster_sizes(p_val_matrix, alpha)
    mask = p_val_matrix < alpha;
    connectedComponents = bwconncomp(mask);
    cluster_sizes = cellfun(@numel, connectedComponents.PixelIdxList);
end

function cluster_sizes = compute_cluster_sizes_z(z_val_matrix, z_thresh)
    mask = z_val_matrix > z_thresh;
    connectedComponents = bwconncomp(mask);
    cluster_sizes = cellfun(@numel, connectedComponents.PixelIdxList);
end

function clusterSizesNull = run_permutations(nPermutations, all_avg_matrices, alpha, use_z, z_thresh)
    clusterSizesNull = cell(nPermutations, 1);  % Preallocate
    for i = 1:nPermutations
        [test_mat, control_mats] = pick_rand_test(all_avg_matrices);
        if use_z
            % Compute z-scored matrix
            z_scored_mat = get_z_mat(test_mat, control_mats);
            clusterSizes = compute_cluster_sizes_z(z_scored_mat, z_thresh);
        else
            % Compute p-value matrix
            p_val_matrix = get_p_mat(test_mat, control_mats);
            clusterSizes = compute_cluster_sizes(p_val_matrix, alpha);
        end
        clusterSizesNull{i} = clusterSizes;
    end
end


function concatenatedVector = concatenate_cluster_sizes(clusterSizesNull)
    concatenatedVector = nan(sum(cellfun(@numel, clusterSizesNull)), 1);
    for i = 1:length(clusterSizesNull)
        concatenatedVector = [concatenatedVector; clusterSizesNull{i}(:)]; % Ensure column vector
    end
end

function sanity_check_matrices(test_mat, control_mats)
    % Check that ctrl_matrix has 3 dimensions
    if ndims(control_mats) ~= 3
        error('ctrl_matrix should have 3 dimensions, but it has %d dimensions.', ndims(control_mats));
    end
    % Check that test_matrix has 2 dimensions
    if ndims(test_mat) ~= 2
        error('test_matrix should have 2 dimensions, but it has %d dimensions.', ndims(test_mat));
    end
    % Check that the first two dimensions of ctrl_matrix match the size of test_matrix
    if size(control_mats, 1) ~= size(test_mat, 1) || size(control_mats, 2) ~= size(test_mat, 2)
        error('The first two dimensions of ctrl_matrix (%d, %d) must match the size of test_matrix (%d, %d).', ...
              size(control_mats, 1), size(control_mats, 2), size(test_mat, 1), size(test_mat, 2));
    end
end

%% single obs plotting
function plot_cluster_size_pdf(real_cluster_sizes, concatenatedVector, nPermutations)
    figure;
    hold on;
    
    % Plot null distribution
    [pdfValues, xValues] = ksdensity(concatenatedVector, 'Bandwidth', 1);
    plot(xValues, pdfValues, 'LineWidth', 2, 'Color', 'Blue', 'DisplayName', ['Permuted Null Distribution n=' num2str(nPermutations)]);
    
    % Plot real data
    [pdfValues, xValues] = ksdensity(real_cluster_sizes(:), 'Bandwidth', 1);
    plot(xValues, pdfValues, 'LineWidth', 2.5, 'Color', 'Black', 'DisplayName', 'Actual');
    
    % Labels and legend
    xlabel('Cluster Size');
    ylabel('Probability Density');
    title('Cluster Size PDF');
    grid on;
    legend('show', 'Location', 'NorthEast');
    
    hold off;
end

function visualize_significant_clusters(real_p_val_matrix, real_cluster_sizes, threshold, alpha)
    figure;
    imagesc(real_p_val_matrix);
    hold on;
    colormap('gray');
    h = colorbar;
    ylabel(h, 'P-Value');
    caxis([0 alpha]);

    % Find connected components
    real_mask = real_p_val_matrix < alpha;
    connectedComponents = bwconncomp(real_mask);
    significant_clusters = real_cluster_sizes > threshold;
    
    % Plot clusters
    for i = 1:length(real_cluster_sizes)
        cluster_indices = connectedComponents.PixelIdxList{i};
        [row, col] = ind2sub(size(real_p_val_matrix), cluster_indices);
        
        if significant_clusters(i)
            plot(col, row, 'r.', 'MarkerSize', 10);
        else
            plot(col, row, 'b.', 'MarkerSize', 5);
        end
    end
    
    title('Clusters with Significant Size');
    xlabel('Maintenance Indices');
    ylabel('Encoding Indices');
    legend('Significant Cluster', 'Non-Significant Cluster', 'Location', 'Best');
    grid on;
    hold off;
end

%% old functions
% function for choosing a random matrix as the test
function [test_mat, control_mats] = pick_rand_test(all_avg_matrices)
    % Step 1: Generate a random index to select the test matrix
    random_index = randi(size(all_avg_matrices, 3));
    % Step 2: Extract the randomly selected matrix as test_mat
    test_mat = all_avg_matrices(:, :, random_index);
    % Step 3: Remove the selected test matrix from the original matrix to get control_mats
    control_mats = all_avg_matrices(:, :, [1:random_index-1, random_index+1:end]);
end

function [z_values, p_values] = get_z_p_mat(test_mat, control_mats)
    % Pre-allocate outputs
    z_values = zeros(size(test_mat));
    p_values = zeros(size(test_mat));
    
    % Compute mean and std along the third dimension in one step
    dist_mean = mean(control_mats, 3);
    dist_std = std(control_mats, 0, 3); % 0 specifies normalization by N-1
    
    % Vectorized z-score computation
    z_values = (test_mat - dist_mean) ./ dist_std;
    
    % Reshape test_mat to broadcast over the third dimension
    test_mat_expanded = reshape(test_mat, size(test_mat, 1), size(test_mat, 2), 1);
    
    % Vectorized percentile computation
    p_values = sum(control_mats > test_mat_expanded, 3) ./ size(control_mats, 3);
end

% legacy ref for above
% function [z_values, p_values] = get_z_p_mat(test_mat, control_mats) % can break into 2 separate functions
%     z_values = zeros(size(test_mat)); % Initialize z_values with the same size as test_matrix
%     p_values = zeros(size(test_mat)); % Initialize p_values with the same size as test_matrix
%     for ide = 1:size(control_mats, 1) % encoding indices
%         for idm = 1:size(control_mats, 2) % maintenance indices
%                 distribution = squeeze(control_mats(ide, idm, :)); % Extract distribution across the third dimension
%                 test_value = test_mat(ide, idm); % Corresponding test matrix value
%                 % Calculate mean and standard deviation of the control distribution
%                 dist_mean = mean(distribution);
%                 dist_std = std(distribution);
%                 % Compute z-score
%                 z_values(ide, idm) = (test_value - dist_mean) / dist_std;
%                 % compute p-value
%                 percentile = sum(distribution > test_value) / length(distribution); % Percent of control values greater than the test value
%                 p_values(ide, idm) = percentile;
%         end
%     end
% end

% different compartmentalization for above
% % function for z-scoring the test by the control distribution
% function p_values = get_p_mat(test_mat, control_mats)
%     % Pre-allocate p_values
%     p_values = zeros(size(test_mat));
% 
%     % Reshape test_mat to broadcast over the third dimension
%     test_mat_expanded = reshape(test_mat, size(test_mat, 1), size(test_mat, 2), 1);
% 
%     % Vectorized percentile computation
%     p_values = sum(control_mats > test_mat_expanded, 3) ./ size(control_mats, 3);
% end
% 
% 
% function z_values = get_z_mat(test_mat, control_mats)
%     % Pre-allocate z_values
%     z_values = zeros(size(test_mat));
% 
%     % Compute mean and std along the third dimension in one step
%     dist_mean = mean(control_mats, 3);
%     dist_std = std(control_mats, 0, 3); % 0 specifies normalization by N-1
% 
%     % Vectorized z-score computation
%     z_values = (test_mat - dist_mean) ./ dist_std;
% end

%% by anat plotting

function clusters_by_uniq_anat = sort_by_uniq_anat(clusters_by_anat)
    % Extract unique anatomical labels
    unique_anat = unique(clusters_by_anat.anat);
    n_anat = length(unique_anat);

    % Initialize structure
    clusters_by_uniq_anat = struct();
    clusters_by_uniq_anat.anat = unique_anat;
    clusters_by_uniq_anat.num_sig_test_clusters = cell(n_anat, 1);
    clusters_by_uniq_anat.sig_test_cluster_sizes = cell(n_anat, 1);
    clusters_by_uniq_anat.real_test_cluster_sizes = cell(n_anat, 1);

    % Aggregate clusters by unique anatomy
    for anat_idx = 1:n_anat
        anat_label = unique_anat{anat_idx};
        % Find indices matching the current anatomical label
        matching_indices = strcmp(clusters_by_anat.anat, anat_label);

        % Aggregate data for this anatomical label
        clusters_by_uniq_anat.num_sig_test_clusters{anat_idx} = ...
            vertcat(clusters_by_anat.num_sig_test_clusters{matching_indices});
        clusters_by_uniq_anat.sig_test_cluster_sizes{anat_idx} = ...
            vertcat(clusters_by_anat.sig_test_cluster_sizes{matching_indices});
        clusters_by_uniq_anat.real_test_cluster_sizes{anat_idx} = ...
            vertcat(clusters_by_anat.real_test_cluster_sizes{matching_indices});
    end
end

% legacy reference for above
% function [clusterSizes_ByAnatLabel] = get_ClusterSizes_byAnatLabel(anatLabel_by_filePairs, clusterSizes_by_all_file_pairs)
%     % Find unique anatomical labels
%     uniqueAnatLabels = unique(anatLabel_by_filePairs);
% 
%     % Initialize a cell array to hold aggregated cluster sizes per anatomy
%     clusterSizes_ByAnatLabel = cell(length(uniqueAnatLabels), 1);
% 
%     % Loop through each file pair and aggregate cluster sizes by anatomical label
%     for i = 1:length(clusterSizes_by_all_file_pairs)
%         anatIndex = find(strcmp(anatLabel_by_filePairs{i}, uniqueAnatLabels));
% 
%         % Ensure clusterSizes is a column vector before concatenation
%         if isrow(clusterSizes_by_all_file_pairs{i})
%             clusterSizes_by_all_file_pairs{i} = clusterSizes_by_all_file_pairs{i}';
%         end
% 
%         % Concatenate cluster sizes for the current anatomical label
%         clusterSizes_ByAnatLabel{anatIndex} = [clusterSizes_ByAnatLabel{anatIndex}; clusterSizes_by_all_file_pairs{i}];
%     end
% end

function plot_cluster_sizes_by_uniq_anat(clusters_by_uniq_anat)
    % Extract unique anatomical labels and corresponding cluster sizes
    unique_anat = clusters_by_uniq_anat.anat;
    n_anat = length(unique_anat);
    cluster_sizes = clusters_by_uniq_anat.real_test_cluster_sizes;

    % Prepare data for boxplot
    all_cluster_sizes = [];
    group_labels = [];
    for anat_idx = 1:n_anat
        all_cluster_sizes = [all_cluster_sizes; cluster_sizes{anat_idx}];
        group_labels = [group_labels; repmat(unique_anat(anat_idx), length(cluster_sizes{anat_idx}), 1)];
    end

    % Create the boxplot
    figure;
    boxplot(all_cluster_sizes, group_labels, 'LabelOrientation', 'inline');
    xlabel('Anatomical Labels');
    ylabel('Cluster Size');
    title('Cluster Sizes Grouped by Anatomical Regions');
    set(gca, 'XTickLabelRotation', 45); % Rotate labels for better visibility
    grid on;
end
% legacy ref for above
% function plotClusterSizes_byAnatLabel(clusterSizes_ByAnatLabel, totalItems_byAnat, uniqueAnatLabels, titleprefix)
%     % Prepare data for boxplot
%     allClusterSizes = cat(1, clusterSizes_ByAnatLabel{:}); % Concatenate all cluster sizes
%     groupLabels = repelem(uniqueAnatLabels, cellfun(@length, clusterSizes_ByAnatLabel)); % Create group labels for boxplot
% 
%     % Create new labels that include the total items (n=)
%     updatedAnatLabels = strcat(uniqueAnatLabels, ' (n=', string(totalItems_byAnat), ')');
% 
%     % Create the box plot
%     figure;
%     boxplot(allClusterSizes, groupLabels, 'LabelOrientation', 'inline');
% 
%     % Update x-axis labels with total items
%     set(gca, 'XTickLabel', updatedAnatLabels);
% 
%     xlabel('Anatomy');
%     ylabel('Cluster Size');
%     title([titleprefix ' Cluster Size Distributions Across Anatomies']);
% end

function plot_percent_sig_by_uniq_anat(clusters_by_uniq_anat)
error("Not Implemented")
end
% legacy ref
% function plotPercentageSignificant_byAnatLabel(percentSignificant_ByAnatLabel, totalItems_byAnat, uniqueAnatLabels)
%     % Sort the data by the magnitude of percentSignificant_ByAnatLabel
%     [sortedPercent, sortIdx] = sort(percentSignificant_ByAnatLabel, 'descend');
%     sortedAnatLabels = uniqueAnatLabels(sortIdx);
%     sortedAnatLabels = strrep(sortedAnatLabels, '_', ' ');
%     sortedTotalItems = totalItems_byAnat(sortIdx);
% 
%     % Plotting
%     figure('Position', [100, 100, 1600, 600]); % Increase figure width for more space
%     bar(sortedPercent);
% 
%     % Set X-ticks and ensure all labels are shown
%     xticks(1:length(sortedAnatLabels));
%     set(gca, 'XTickLabel', sortedAnatLabels, 'XTickLabelRotation', 80, 'FontSize', 8);
% 
%     xlabel('Anatomical Labels');
%     ylabel('Percentage of Ch/Item combinations with >0 Significant Clusters');
%     title('Significant Replay related to item performance by Region (Z-scored mean itemCorrect EMS to itemWrong EMS.)');
%     grid on;
%     ylim([0 100]);
% 
%     % Display total items per label on the plot (as text above bars)
%     xtips1 = 1:length(sortedPercent); 
%     ytips1 = sortedPercent;
%     labels1 = strcat('n=', string(sortedTotalItems));
%     text(xtips1, ytips1, labels1, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8);
% end

function plot_mean_sig_cluster_size_by_anat(clusters_by_uniq_anat)
error("not implemented")
end

function plot_max_sig_cluster_size_by_anat(clusters_by_uniq_anat)
error("not implemented")
end

function plot_num_sig_clusters_by_anat(clusters_by_uniq_anat)
error("not implemented")
end