% if temporally generalized use this witohut parfor
function [t_values, p_values] = calc_ttest(diff)
    % tic;
    dim1 = size(diff,1);
    dim2 = size(diff,2);

    t_values = zeros(dim1, dim2);
    p_values = zeros(dim1, dim2);

    if dim1 > 1 % not temporally generalized size=(E,M,K) use parfor.
        parfor i = 1:dim1
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
        % fprintf("finished parfor t-test in %s\n", num2str(toc))

    else % temporally generalized size=(1,1,K)
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
        % fprintf("finished for-loop t-test in %s\n", num2str(toc))
    end

    if any(isnan(p_values), "all") || any(isnan(t_values), "all")
        error("t_values or p_values from t-test contains nan")
    end
end

% gpus overhead is slower it seems
% function [t_values, p_values] = calc_ttest(diff)
%     tic;
%     % diff = gpuArray(diff);
%     [dim1, dim2, ~] = size(diff); % Get the dimensions
%     t_values = gpuArray.zeros(dim1, dim2);  % Initialize t-values
%     p_values = gpuArray.zeros(dim1, dim2);  % Initialize p-values
% 
%     % Precompute the list of indices in 1D (linear indices)
%     indices = gpuArray(1:(dim1 * dim2));
% 
%     % Flatten the diff matrix to avoid ind2sub
%     diff_flat = gpuArray(reshape(diff, [], size(diff, 3)));  % Flatten to [dim1*dim2, n_samples] and move to GPU
% 
% 
%     % Use parfor to iterate over the linear index list
%     parfor k = indices
%         % data_diff = gpuArray(diff_flat(k, :));  % This directly accesses the data for (i, j)
% 
%         % [i, j] = ind2sub([dim1, dim2], k);  % Convert linear index to (i, j)
% 
%         % Extract the distribution for the current time-time pair
%         % data_diff = squeeze(diff(i, j, :));
% 
%         % Perform one-sample t-test against 0
%         [~, p, ~, stats] = ttest(diff_flat(k, :));
% 
%         % Store the t-value and p-value in linear indices
%         t_values(k) = stats.tstat;
%         p_values(k) = p;
%     end
% 
%     % Convert the 1D results back to 2D
%     t_values = reshape(t_values, dim1, dim2);
%     p_values = reshape(p_values, dim1, dim2);
% 
%     % Optionally, gather results back to CPU
%     t_values = gather(t_values);
%     p_values = gather(p_values);
%     fprintf("finished gpu t-test in %s\n", num2str(toc))
% end
