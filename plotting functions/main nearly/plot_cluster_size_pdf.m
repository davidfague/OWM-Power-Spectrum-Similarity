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