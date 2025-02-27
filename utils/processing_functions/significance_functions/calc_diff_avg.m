function diff = calc_diff_avg(WI, BI, num_means)
    if size(WI, 1) > 1
        % tic;
    
        % Move input matrices to GPU
        try
            WI = gpuArray(WI);
            BI = gpuArray(BI);
            has_gpus = true;
        catch
            has_gpus = false;
        end
    
        % Get dimensions
        [t1_size, t2_size, num_WI] = size(WI);
        num_BI = size(BI, 3);
        min_size = min(num_WI, num_BI);
    
        % Precompute random indices on the GPU
        if num_WI < num_BI
            if has_gpus
                random_indices_WI = gpuArray(repmat(1:num_WI, [num_means, 1]));
                random_indices_BI = gpuArray(randi(num_BI, num_means, num_WI));
            else
                random_indices_WI = repmat(1:num_WI, [num_means, 1]);
                random_indices_BI = randi(num_BI, num_means, num_WI);
            end
        else
            if has_gpus
                random_indices_WI = gpuArray(randi(num_WI, num_means, num_BI));
                random_indices_BI = gpuArray(repmat(1:num_BI, [num_means, 1]));
            else
                random_indices_WI = randi(num_WI, num_means, num_BI);
                random_indices_BI = repmat(1:num_BI, [num_means, 1]);
            end
        end
    
        % Flatten WI and BI matrices for processing
        WI = reshape(WI, [], num_WI); % (t1 * t2) x num_WI
        BI = reshape(BI, [], num_BI); % (t1 * t2) x num_BI
    
        % Preallocate result matrices on GPU
        if has_gpus
            WI_means = gpuArray.zeros(t1_size * t2_size, num_means, 'single');
            BI_means = gpuArray.zeros(t1_size * t2_size, num_means, 'single');
        else
            WI_means = zeros(t1_size * t2_size, num_means, 'single');
            BI_means = zeros(t1_size * t2_size, num_means, 'single');
        end
    
        % % Compute means across random indices in batches
        for k = 1:num_means
            WI_means(:, k) = mean(WI(:, random_indices_WI(k, :)), 2);
            BI_means(:, k) = mean(BI(:, random_indices_BI(k, :)), 2);
        end
    
        % Compute difference
        diff = reshape(WI_means - BI_means, t1_size, t2_size, num_means);
    
        % Bring results back to CPU
        % diff = gather(diff);
    
        % fprintf("Finished in %.2f seconds\n", toc);
    else
        % use this if temporally generalized
            % Move the input matrices to GPU % in this case overhead is not worth since 38 x 380 time dimension becomes 1x1
            % WI = gpuArray(WI);
            % BI = gpuArray(BI);
        
            % Dimensions of the input matrices
            num_WI = size(WI, 3);
            num_BI = size(BI, 3);
        
            % Flatten WI and BI along their 3rd dimension (e.g., reshape them into 2D)
            WI_flat = reshape(WI, [], num_WI);  % Flatten WI to (time*time) x num_WI
            BI_flat = reshape(BI, [], num_BI);  % Flatten BI to (time*time) x num_BI
        
            % Calculate the minimum size for the 3rd dimension
            min_size = min(num_WI, num_BI);
        
            % Generate random sets of indices to select from the 3rd dimension
            if num_WI < num_BI
                % Generate `num_pairs` sets of random indices for WI (1:num_WI) and BI (min_size random from 1:num_BI)
                random_indices_WI = arrayfun(@(x) 1:num_WI, 1:num_means, 'UniformOutput', false);
                random_indices_BI = arrayfun(@(x) randperm(num_BI, min_size), 1:num_means, 'UniformOutput', false);
            else
                % Generate `num_pairs` sets of random indices for WI (min_size random from 1:num_WI) and BI (1:num_BI)
                random_indices_WI = arrayfun(@(x) randperm(num_WI, min_size), 1:num_means, 'UniformOutput', false);
                random_indices_BI = arrayfun(@(x) 1:num_BI, 1:num_means, 'UniformOutput', false);
            end
        
            % this method chooses random, repeating
            % random_indices_WI = randi(num_WI, [num_pairs, min_size]);  % num_pairs sets of random indices for WI
            % random_indices_BI = randi(num_BI, [num_pairs, min_size]);  % num_pairs sets of random indices for BI
        
            % Preallocate the output difference matrix
            diff = zeros(size(WI, 1), size(WI, 2), num_means);
        
            % Calculate differences for the selected combinations
        
            for k = 1:num_means
                % Select the random indices for WI and BI
                % indices_WI = random_indices_WI(k, :); % for randi, repeating
                % indices_BI = random_indices_BI(k, :);
        
                indices_WI = random_indices_WI{k}; % for randperm
                indices_BI = random_indices_BI{k};
        
                % Extract the data from the flattened matrices
                WI_data = WI_flat(:, indices_WI);  % Select the data for WI from the flattened matrix
                BI_data = BI_flat(:, indices_BI);  % Select the data for BI from the flattened matrix
        
                % % Compute mean across the selected indices for WI and BI
                % WI_mean = mean(WI(:, :, indices_WI), 3);
                % BI_mean = mean(BI(:, :, indices_BI), 3);
        
                % Compute mean across the selected indices for WI and BI
                WI_mean = mean(WI_data, 2);  % Mean across the selected indices (2nd dimension)
                BI_mean = mean(BI_data, 2);  % Mean across the selected indices (2nd dimension)
        
                % Compute the difference
                % diff(:, :, k) = WI_mean - BI_mean;
                diff(:, :, k) = reshape(WI_mean - BI_mean, size(WI, 1), size(WI, 2));
        
            end
            % diff = gather(diff);
            % fprintf("finished gpu diffs in %s\n", num2str(toc))
    end
end


% trying to use large matrices on gpu instead of loops
% function diff = calc_diff_avg(WI, BI, num_means)
%     tic
%     % Inputs:
%     % WI - The within-item matrix (size: time x time x num_WI)
%     % BI - The between-item matrix (size: time x time x num_BI)
%     % num_means - number of means to compute (each mean is across min_size
%     % combinations of WI and BI then the first WI and first BI mean are
%     % subtracted, seconds, thirds, etc.
% 
%     % Move the input matrices to GPU
%     WI = gpuArray(WI);
%     BI = gpuArray(BI);
% 
%     % Dimensions of the input matrices
%     num_WI = size(WI, 3);
%     num_BI = size(BI, 3);
%     t1_size = size(WI, 1);
%     t2_size = size(WI, 2);
% 
%     % Flatten time dimensions of WI and BI along their 3rd dimension (e.g.,
%     % reshape them into 2D) (times, combinations) (combined time dimensions)
%     WI = reshape(WI, [], num_WI);  % Flatten WI to (time*time) x num_WI
%     BI = reshape(BI, [], num_BI);  % Flatten BI to (time*time) x num_BI
% 
%     % Calculate the minimum size for the 3rd dimension
%     min_size = min(num_WI, num_BI);
% 
%     % Generate random sets of indices for WI and BI (as numeric arrays)
%     % (these slice the combinations dimension 1000 times) (we have to
%     % downsample one of the datasets so they have the same number of
%     % combinations
%     if num_WI < num_BI
%         % Generate `num_means` sets of random indices for WI (1:num_WI) and BI (min_size random from 1:num_BI)
%         random_indices_WI = gpuArray(repmat(1:num_WI,[num_means, 1]));
%         WI = mean(reshape(WI(:,random_indices_WI(:)), size(WI, 1), size(random_indices_WI, 1), []),3);
%         clear random_indices_WI
% 
%         % random_indices_BI = gpuArray(randi(num_BI, [num_pairs, num_WI]));
%         random_indices_BI = arrayfun(@(x) randperm(num_BI, num_WI), 1:num_means, 'UniformOutput', false);
%         random_indices_BI = gpuArray(cell2mat(random_indices_BI'));
%         BI = mean(reshape(BI(:,random_indices_BI(:)), size(BI, 1), size(random_indices_BI, 1), []),3);
%         clear random_indices_BI
%     else
%         % Generate `num_means` sets of random indices for WI (min_size random from 1:num_WI) and BI (1:num_BI)
%         random_indices_WI = arrayfun(@(x) randperm(num_WI, num_BI), 1:num_means, 'UniformOutput', false);
%         random_indices_WI = gpuArray(cell2mat(random_indices_WI'));
%         WI = reshape(WI(:,random_indices_WI(:)), size(WI, 1), size(random_indices_WI, 1), []);
% 
%         random_indices_BI = gpuArray(repmat(1:num_BI,[num_means, 1]));  % Create a numeric array directly  % Directly assign a range for BI
%         BI = reshape(BI(:,random_indices_BI(:)), size(BI, 1), size(random_indices_BI, 1), []);
%         clear random_indices_BI
%     end
% 
%     % Compute the means across the third dimension (num_pairs)
%     % WI_mean = mean(WI_combined, 2);  % Compute mean across pairs (time*time x num_averages)
%     % BI_mean = mean(BI_combined, 2);  % Compute mean across pairs (time*time x num_averages)
% 
%     % Reshape to match original dimensions and compute the difference
%     diff = reshape(WI - BI, t1_size, t2_size, num_means);
% 
%     diff = gather(diff);  % Bring the result back to CPU
%     fprintf("finished gpu diffs in %s\n", num2str(toc))
% end

% bad
% function diff = calc_diff_avg(WI, BI, num_means)
%     tic
%     % Inputs:
%     % WI - The within-item matrix (size: time x time x num_WI)
%     % BI - The between-item matrix (size: time x time x num_BI)
%     % num_means - number of means to compute (each mean is across min_size
%     % combinations of WI and BI then the first WI and first BI mean are
%     % subtracted, seconds, thirds, etc.
% 
%     % Move the input matrices to GPU
%     WI = gpuArray(WI);
%     BI = gpuArray(BI);
% 
%     % Dimensions of the input matrices
%     num_WI = size(WI, 3);
%     num_BI = size(BI, 3);
%     t1_size = size(WI, 1);
%     t2_size = size(WI, 2);
% 
%     % Flatten time dimensions of WI and BI along their 3rd dimension (e.g.,
%     % reshape them into 2D) (times, combinations) (combined time dimensions)
%     WI_flat = reshape(WI, [], num_WI);  % Flatten WI to (time*time) x num_WI
%     BI_flat = reshape(BI, [], num_BI);  % Flatten BI to (time*time) x num_BI
%     clear BI WI
% 
%     % Calculate the minimum size for the 3rd dimension
%     min_size = min(num_WI, num_BI);
% 
%     % Generate random sets of indices for WI and BI (as numeric arrays)
%     % (these slice the combinations dimension 1000 times) (we have to
%     % downsample one of the datasets so they have the same number of
%     % combinations
%     if num_WI < num_BI
%         % Generate `num_means` sets of random indices for WI (1:num_WI) and BI (min_size random from 1:num_BI)
%         random_indices_WI = gpuArray(repmat(1:num_WI,[num_means, 1]));  % Create a numeric array directly
%         WI_combined = mean(reshape(WI_flat(:,random_indices_WI(:)), size(WI_flat, 1), size(random_indices_WI, 1), []),3);
%         clear WI_flat random_indices_WI
% 
%         % random_indices_BI = gpuArray(randi(num_BI, [num_pairs, num_WI]));  % Create a numeric array directly
%         random_indices_BI = arrayfun(@(x) randperm(num_BI, num_WI), 1:num_means, 'UniformOutput', false);
%         random_indices_BI = gpuArray(cell2mat(random_indices_BI'));
%         BI_combined = mean(reshape(BI_flat(:,random_indices_BI(:)), size(BI_flat, 1), size(random_indices_BI, 1), []),3);
%         clear BI_flat random_indices_BI
%     else
%         % Generate `num_means` sets of random indices for WI (min_size random from 1:num_WI) and BI (1:num_BI)
%         random_indices_WI = arrayfun(@(x) randperm(num_WI, num_BI), 1:num_means, 'UniformOutput', false);
%         random_indices_WI = gpuArray(cell2mat(random_indices_WI'));
%         WI_combined = reshape(WI_flat(:,random_indices_WI(:)), size(WI_flat, 1), size(random_indices_WI, 1), []);
%         clear WI_flat random_indices_WI
% 
%         random_indices_BI = gpuArray(repmat(1:num_BI,[num_means, 1]));  % Create a numeric array directly  % Directly assign a range for BI
%         BI_combined = reshape(BI_flat(:,random_indices_BI(:)), size(BI_flat, 1), size(random_indices_BI, 1), []);
%         clear BI_flat random_indices_BI
%     end
% 
%     % Compute the means across the third dimension (num_pairs)
%     % WI_mean = mean(WI_combined, 2);  % Compute mean across pairs (time*time x num_averages)
%     % BI_mean = mean(BI_combined, 2);  % Compute mean across pairs (time*time x num_averages)
% 
%     % Reshape to match original dimensions and compute the difference
%     diff = reshape(WI_combined - BI_combined, t1_size, t2_size, num_means);
% 
%     diff = gather(diff);  % Bring the result back to CPU
%     fprintf("finished gpu diffs in %s\n", num2str(toc))
% end
% 
% function diff = calc_diff_avg(WI, BI, num_pairs)
%     % tic
%     % Inputs:
%     % WI - The within-item matrix (size: time x time x num_WI)
%     % BI - The between-item matrix (size: time x time x num_BI)
%     % num_pairs - number of pairs to compute (randomly selected if less than total_pairs)
% 
%     % Move the input matrices to GPU
%     WI = gpuArray(WI);
%     BI = gpuArray(BI);
% 
%     % Dimensions of the input matrices
%     num_WI = size(WI, 3);
%     num_BI = size(BI, 3);
% 
%     % Flatten WI and BI along their 3rd dimension (e.g., reshape them into 2D)
%     WI_flat = reshape(WI, [], num_WI);  % Flatten WI to (time*time) x num_WI
%     BI_flat = reshape(BI, [], num_BI);  % Flatten BI to (time*time) x num_BI
% 
%     % Calculate the minimum size for the 3rd dimension
%     min_size = min(num_WI, num_BI);
% 
%     % % Generate random sets of indices to select from the 3rd dimension
%     % if num_WI < num_BI
%     %     % Generate `num_pairs` sets of random indices for WI (1:num_WI) and BI (min_size random from 1:num_BI)
%     %     random_indices_WI = arrayfun(@(x) 1:num_WI, 1:num_pairs, 'UniformOutput', false);
%     %     random_indices_BI = arrayfun(@(x) randperm(num_BI, min_size), 1:num_pairs, 'UniformOutput', false);
%     % else
%     %     % Generate `num_pairs` sets of random indices for WI (min_size random from 1:num_WI) and BI (1:num_BI)
%     %     random_indices_WI = arrayfun(@(x) randperm(num_WI, min_size), 1:num_pairs, 'UniformOutput', false);
%     %     random_indices_BI = arrayfun(@(x) 1:num_BI, 1:num_pairs, 'UniformOutput', false);
%     % end
% 
%     % Generate random sets of indices for WI and BI (as numeric arrays)
%     if num_WI < num_BI
%         % Generate `num_pairs` sets of random indices for WI (1:num_WI) and BI (min_size random from 1:num_BI)
%         random_indices_WI = gpuArray(randi(num_WI, [num_pairs, 1]));  % Create a numeric array directly
%         random_indices_BI = gpuArray(randi(num_BI, [num_pairs, min_size]));  % Create a numeric array directly
%     else
%         % Generate `num_pairs` sets of random indices for WI (min_size random from 1:num_WI) and BI (1:num_BI)
%         random_indices_WI = gpuArray(randi(num_WI, [num_pairs, min_size]));  % Create a numeric array directly
%         random_indices_BI = gpuArray(1:num_BI);  % Directly assign a range for BI
%     end
% 
%     % Precompute slices for each k (each pair) to avoid broadcasting in parfor
%     WI_data_precomputed = cell(1, num_pairs);
%     BI_data_precomputed = cell(1, num_pairs);
% 
%     for k = 1:num_pairs
%         % Select and precompute the data slices for each pair (WI and BI)
%         indices_WI = random_indices_WI(k);
%         indices_BI = random_indices_BI(k);
% 
%         % Precompute the data for the current pair
%         WI_data_precomputed{k} = WI_flat(:, indices_WI);  % Select the data for WI from the flattened matrix
%         BI_data_precomputed{k} = BI_flat(:, indices_BI);  % Select the data for BI from the flattened matrix
%     end
% 
%     % Preallocate the output difference matrix
%     diff = gpuArray.zeros(size(WI, 1), size(WI, 2), num_pairs);
% 
%     % Calculate differences for the selected combinations in parfor
%     parfor k = 1:num_pairs
%         % Access precomputed data for the current iteration
%         % WI_data = WI_data_precomputed{k};
%         % BI_data = BI_data_precomputed{k};
% 
%         % Compute mean across the selected indices for WI and BI
%         WI_mean = mean(WI_data_precomputed{k}, 2);  % Mean across the selected indices (2nd dimension)
%         BI_mean = mean(BI_data_precomputed{k}, 2);  % Mean across the selected indices (2nd dimension)
% 
%         % Reshape to match original dimensions and compute the difference
%         diff(:, :, k) = reshape(WI_mean - BI_mean, size(WI, 1), size(WI, 2));
%     end
% 
%     diff = gather(diff);  % Bring the result back to CPU
%     % fprintf("finished gpu diffs in %s\n", num2str(toc))
% end


% % use this if temporally generalized
% function diff = calc_diff_avg(WI, BI, num_pairs)
%     % tic
%     % Inputs:
%     % WI - The within-item matrix (size: time x time x num_WI)
%     % BI - The between-item matrix (size: time x time x num_BI)
%     % num_pairs - number of pairs to compute (randomly selected if less than total_pairs)
% 
%     % Move the input matrices to GPU
%     % WI = gpuArray(WI);
%     % BI = gpuArray(BI);
% 
%     % Dimensions of the input matrices
%     num_WI = size(WI, 3);
%     num_BI = size(BI, 3);
% 
%     % Flatten WI and BI along their 3rd dimension (e.g., reshape them into 2D)
%     WI_flat = reshape(WI, [], num_WI);  % Flatten WI to (time*time) x num_WI
%     BI_flat = reshape(BI, [], num_BI);  % Flatten BI to (time*time) x num_BI
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
% 
%     for k = 1:num_pairs
%         % Select the random indices for WI and BI
%         % indices_WI = random_indices_WI(k, :); % for randi, repeating
%         % indices_BI = random_indices_BI(k, :);
% 
%         indices_WI = random_indices_WI{k}; % for randperm
%         indices_BI = random_indices_BI{k};
% 
%         % Extract the data from the flattened matrices
%         WI_data = WI_flat(:, indices_WI);  % Select the data for WI from the flattened matrix
%         BI_data = BI_flat(:, indices_BI);  % Select the data for BI from the flattened matrix
% 
%         % % Compute mean across the selected indices for WI and BI
%         % WI_mean = mean(WI(:, :, indices_WI), 3);
%         % BI_mean = mean(BI(:, :, indices_BI), 3);
% 
%         % Compute mean across the selected indices for WI and BI
%         WI_mean = mean(WI_data, 2);  % Mean across the selected indices (2nd dimension)
%         BI_mean = mean(BI_data, 2);  % Mean across the selected indices (2nd dimension)
% 
%         % Compute the difference
%         % diff(:, :, k) = WI_mean - BI_mean;
%         diff(:, :, k) = reshape(WI_mean - BI_mean, size(WI, 1), size(WI, 2));
% 
%     end
%     % diff = gather(diff);
%     % fprintf("finished gpu diffs in %s\n", num2str(toc))
% end