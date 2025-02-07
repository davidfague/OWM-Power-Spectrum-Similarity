function fig = plot_t_values(real_value, null_values)
    % Plot the histogram of null values
    fig = figure;
    histogram(null_values, 30, 'Normalization', 'pdf');  % 30 bins, normalized to probability density
    hold on;
    
    % Add a vertical line at the real value
    xline(real_value, 'r-', 'LineWidth', 2);  % Red line for the real value
    
    % Add a label for the real value line with p-value calculation
    p_value = mean(real_value < null_values);  % p-value as the proportion of null values less than the real value
    text(real_value, 0.1, sprintf('p = %.3f', p_value), 'Color', 'red', 'FontSize', 12);
    
    % Add titles and labels
    title('Distribution of Null Values with Real Value Line');
    xlabel('Value');
    ylabel('Probability Density');
    grid on;
    hold off;
end