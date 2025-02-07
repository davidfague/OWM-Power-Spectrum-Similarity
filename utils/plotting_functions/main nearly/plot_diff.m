% function plot_diff(diff)
% figure;
% histogram(diff, 'Normalization', 'pdf');
% title('PDF of diff');
% xlabel('Value');
% ylabel('Probability Density');
% end

function fig = plot_diff(diff1, diff2, final_p, plot_params)
% note that real diff should come first
    fig = figure;
    hold on;  % Hold the plot to overlay both histograms
    legend('AutoUpdate','on');
    histogram(diff1, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'DisplayName', 'Real Differences');  % Plot the first diff with some transparency
    histogram(diff2, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'DisplayName', 'Null Differences');  % Plot the second diff with some transparency
    hold off;  % Release the plot hold
    
    if plot_params.average_diff
        title(sprintf('PDF of WI-BI (average trials(N=1000)) (p < %s)', num2str(abs(1-final_p))));
    else
        title(sprintf('PDF of WI-BI (individual trials(N=10000)) (p < %s)', num2str(abs(1-final_p))));
    end
    xlabel('Value');
    ylabel('Probability Density');
    % legend({'diff1', 'diff2'});  % Add a legend to distinguish the two PDFs
end
