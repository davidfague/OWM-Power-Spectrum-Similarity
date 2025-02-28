function p_values = assess_significance(real_sums, surrogate_sums, cluster_type)
    % Compute p-values by comparing real cluster sums to surrogate distribution
    p_values = zeros(size(real_sums));
    if strcmp(cluster_type, 'positive')
        % Positive clusters: Compare to max surrogate sums
        for i = 1:length(real_sums)
            p_values(i) = mean(surrogate_sums >= real_sums(i)); % right tail
        end
    elseif strcmp(cluster_type, 'negative')
        % Negative clusters: Compare to min surrogate sums
        for i = 1:length(real_sums)
            p_values(i) = mean(surrogate_sums <= real_sums(i)); % left tail
        end
    else
        error('Invalid cluster type. Use "positive" or "negative".');
    end
end