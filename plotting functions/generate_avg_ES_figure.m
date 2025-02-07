function fig = generate_avg_ES_figure(avg_ES, folder_to_analyze)
    % can turn this into function and clean up
    fig = figure('WindowState','maximized');%, 'Visible', 'off');%, 'Visible', 'off');
    title_str = strrep(folder_to_analyze, '_', ' ');
    title_str = strcat('Average Within-Trial Whole-Trial Encoding Similarity', title_str);
    title(title_str, 'Interpreter', 'none')
    imagesc(avg_ES)
    c = colorbar;
    c.Label.String = 'Similarity (rho)';
    clim([0 1])
    title(title_str);
    xlabel('Whole Trial Window ID');
    ylabel('Encoding Window ID');
    % colormap('hot');

    % Plot the vertical lines
    index_adjust = -5; % subtracting 5 actually picks the window id that is centered at the desired time instead of the window id that begins at the desired time.
    xline(25+index_adjust, 'g', 'LineWidth', 2);    % baseline start
    xline(75+index_adjust, 'g', 'LineWidth', 2);    % baseline end
    xline(100+index_adjust, 'b', 'LineWidth', 2);   % fixation end
    xline(150+index_adjust, 'b', 'LineWidth', 2);   % enc1 end
    xline(200+index_adjust, 'b', 'LineWidth', 2);   % enc2 end
    xline(250+index_adjust, 'b', 'LineWidth', 2);   % enc3 end
    xline(650+index_adjust, 'b', 'LineWidth', 2);   % maintenance end
    label_loc = ylim;
    label_loc = label_loc(2);
    label_loc = label_loc + 5;
    % Add text labels at corresponding positions
    text(25+index_adjust, label_loc, 'baseline start', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'g');
    text(75+index_adjust, label_loc, 'baseline end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'g');
    text(100+index_adjust, label_loc, 'fixation end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(150+index_adjust, label_loc, 'enc1 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(200+index_adjust, label_loc, 'enc2 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(250+index_adjust, label_loc, 'enc3 end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');
    text(650+index_adjust, label_loc, 'maint end', 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'b');

end