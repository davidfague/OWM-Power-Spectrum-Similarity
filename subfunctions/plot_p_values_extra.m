function plot_p_values_extra(final_p_by_channels, anat_by_channels)
    figure;
    % Sort by descending order of final_p_by_channels
    [sorted_p, sort_idx] = sort(final_p_by_channels, 'descend');
    sorted_anat = anat_by_channels(sort_idx);

    % Plot the sorted data
    bar(sorted_p)
    xticks(1:length(sorted_anat)); % Ensure all ticks are present
    xticklabels(sorted_anat)
    ylabel('p-value')
    xlabel('Channels')
    title('Sorted p-values by Channel')
    
    % Rotate x-axis labels for better visibility
    xtickangle(90)
end