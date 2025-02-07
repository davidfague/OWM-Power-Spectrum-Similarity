function plot_t_with_p_highlights(t_matrix, p_logical) % for plotting t and overlaying p
    % Validate input matrices
    if ~isequal(size(t_matrix), size(p_logical))
        error('Input matrices must have the same size.');
    end
    
    % Plot the t-values using imagesc
    figure;
    subplot(2,1,1)
    imagesc(t_matrix);
    colormap('jet'); % Use a colormap suitable for t-values
    colorbar; % Add a colorbar for reference
    axis equal tight; % Adjust axis scaling
    title('T-values with Significant Highlights (p < 0.05)');
    xlabel('X-axis');
    ylabel('Y-axis');
    
    subplot(2,1,2)
    imagesc(t_matrix);
    colormap('jet'); % Use a colormap suitable for t-values
    colorbar; % Add a colorbar for reference
    axis equal tight; % Adjust axis scaling
    title('T-values with Significant Highlights (p < 0.05)');
    xlabel('X-axis');
    ylabel('Y-axis');
    % Overlay highlights for significant values (p < 0.05)
    hold on;
    [rows, cols] = find(p_logical); % Get indices of significant points
    scatter(cols, rows, 50, 'k', 'filled', 'MarkerFaceAlpha', 0.6); % Highlight with semi-transparent black dots
    hold off;
end