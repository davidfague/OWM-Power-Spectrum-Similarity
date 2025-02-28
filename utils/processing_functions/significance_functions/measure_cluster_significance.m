function [clusters_info] = measure_cluster_significance(real_p_values, real_t_values, surrogate_p_values, surrogate_t_values)
    alpha = 0.025; % for separating positive & negative clusters % 2 tailed
    significance_threshold = 0.025; % for cluster T sum significance % 2 tailed

    %% gather REAL clusters and cluster's sum(T)
    % threshold the real data, separating positive and negative clusters
    real_positive_mask = (real_p_values < alpha) & (real_t_values > 0);
    real_negative_mask = (real_p_values < alpha) & (real_t_values < 0);

    % identify clusters (connected components)
    real_positive_clusters = bwconncomp(real_positive_mask, 8); % 8-connectivity
    real_negative_clusters = bwconncomp(real_negative_mask, 8);
    % all_real_clusters = [positive_clusters, negative_clusters];

    % calculate cluster sums for real data
    real_positive_cluster_sums = calculate_cluster_sums(real_t_values, real_positive_clusters);
    real_negative_cluster_sums = abs(calculate_cluster_sums(real_t_values, real_negative_clusters));
    % all_real_cluster_sums = [real_positive_cluster_sums, abs(real_negative_cluster_sums)];

    %% gather SURROGATE clusters and cluster's sum(T)
    % preallocate surrogate clusters and cluster sums
    num_surrogates = size(surrogate_t_values, 3);
    surrogate_positive_sums = zeros(num_surrogates, 1);
    surrogate_negative_sums = zeros(num_surrogates, 1);
    % all_surrogate_sums = zeros(num_surrogates, 1);

    surrogate_positive_masks = nan([size(surrogate_p_values, [1 2]), num_surrogates]);
    surrogate_negative_masks = nan([size(surrogate_p_values, [1 2]), num_surrogates]);

    surrogate_positive_clusters = cell(num_surrogates);
    surrogate_negative_clusters = cell(num_surrogates);

    for i = 1:num_surrogates
        % threshold surrogate and separate positive and negative clusters
        surrogate_positive_masks(:,:,i) = (surrogate_p_values(:, :, i) < alpha) & (surrogate_t_values(:, :, i) > 0);
        surrogate_negative_masks(:,:,i) = (surrogate_p_values(:, :, i) < alpha) & (surrogate_t_values(:, :, i) < 0);

        % find clusters
        surrogate_positive_clusters{i} = bwconncomp(surrogate_positive_masks(:,:,i), 8);
        surrogate_negative_clusters{i} = bwconncomp(surrogate_negative_masks(:,:,i), 8);

        % calculate largest positive and negative cluster sum for this surrogate
        if surrogate_positive_clusters{i}.NumObjects == 0
            surrogate_positive_sums(i) = 0;
        else
            surrogate_positive_sums(i) = max(calculate_cluster_sums(surrogate_t_values(:, :, i), surrogate_positive_clusters{i}), [], 'omitnan');
        end

        if surrogate_negative_clusters{i}.NumObjects == 0
            surrogate_negative_sums(i) = 0;
        else
            surrogate_negative_sums(i) = max(abs(calculate_cluster_sums(surrogate_t_values(:, :, i), surrogate_negative_clusters{i})), [], 'omitnan'); % max because doing abs
        end

        % pool only the largest (not correct, should keep separate)
        % if abs(surrogate_positive_sums(i)) > abs(surrogate_negative_sums(i))
        %     all_surrogate_sums(i) = surrogate_positive_sums(i);
        % else
        %     all_surrogate_sums(i) = surrogate_negative_sums(i);
        % end
    end

    %% significance testing

    % compute percent of surrogate largest clusters sum(t)'s less than real
    % sum(t) for each real cluster
    positive_cluster_p_values = assess_significance(real_positive_cluster_sums, surrogate_positive_sums, 'positive');
    negative_cluster_p_values = assess_significance(real_negative_cluster_sums, surrogate_negative_sums, 'positive'); % use positive because taking abs() before.

    % determine significance for clusters
    significant_positive_clusters = positive_cluster_p_values < significance_threshold;
    significant_negative_clusters = negative_cluster_p_values < significance_threshold;
    % significant_clusters = cluster_p_values < significance_threshold;

    % compute the threshold sum for plotting
    positive_sum_threshold = prctile(surrogate_positive_sums, (1-significance_threshold)*100);
    negative_sum_threshold = prctile(surrogate_negative_sums, (1-significance_threshold)*100);
    % negative_sum_threshold = prctile(surrogate_negative_sums, significance_threshold*100);
    % sum_threshold = prctile(abs(all_surrogate_sums), (1-significance_threshold)*100);

    % p_largest_cluster = min([positive_cluster_p_values, negative_cluster_p_values]);
    p_largest_positive_cluster = min(positive_cluster_p_values);
    p_largest_negative_cluster = min(negative_cluster_p_values);
    % final_p = mean(abs(all_surrogate_sums)>max(all_real_cluster_sums));

    % return all variable in a struct
    vars = whos; % Get list of variables in workspace
    clusters_info = struct(); % Initialize an empty struct
    for i = 1:length(vars)
        clusters_info.(vars(i).name) = eval(vars(i).name);
    end
end