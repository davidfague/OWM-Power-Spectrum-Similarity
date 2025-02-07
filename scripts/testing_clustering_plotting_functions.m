%% to confirm these old functions are ones that I want to use

%% Generate synthetic data for testing

% Parameters
num_clusters = 100;       % Number of clusters
cluster_size_range = [1, 50]; % Range of cluster sizes
matrix_size = [100, 100];  % Size of p-value matrix
threshold = 20;           % Cluster size threshold for significance
alpha = 0.05;             % Significance level
nPermutations = 1000;     % Number of permutations for null distribution

% Generate synthetic real cluster sizes
real_cluster_sizes = randi(cluster_size_range, num_clusters, 1);

% Generate synthetic concatenated vector for null distribution
concatenatedVector = randi(cluster_size_range, nPermutations * num_clusters, 1);

% Generate synthetic p-value matrix
real_p_val_matrix = rand(matrix_size);

% Visualize the data
figure;
subplot(1, 2, 1);
histogram(real_cluster_sizes, 'Normalization', 'pdf', 'FaceColor', 'k');
xlabel('Cluster Size');
ylabel('Probability Density');
title('Synthetic Real Cluster Sizes PDF');

subplot(1, 2, 2);
histogram(concatenatedVector, 'Normalization', 'pdf', 'FaceColor', 'b');
xlabel('Cluster Size');
ylabel('Probability Density');
title('Synthetic Null Distribution PDF');

% Test plot_cluster_size_pdf function
plot_cluster_size_pdf(real_cluster_sizes, concatenatedVector, nPermutations);

% Test visualize_significant_clusters function
visualize_significant_clusters(real_p_val_matrix, real_cluster_sizes, threshold, alpha);
