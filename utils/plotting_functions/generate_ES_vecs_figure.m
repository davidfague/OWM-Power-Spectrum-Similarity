
function fig = generate_ES_vecs_figure(ES_file, folder_name)
    fig = figure('WindowState','maximized');%, 'Visible', 'off');%,'Visible', 'off');
    
    title_str = strrep(folder_name, '_', ' ');
    title_str = strcat("Mean Z Power By-Frequency Vectors ", title_str);

    % Calculate the relative widths based on the data sizes
    encoding_data_size = size(ES_file.avg_window1_mean_PS_vectors, 1); % Encoding x-axis size (40)
    trial_data_size = size(ES_file.avg_window2_mean_PS_vectors, 1); % Whole trial x-axis size (891)
    total_size = encoding_data_size + trial_data_size;
 
    encoding_width = encoding_data_size / total_size; % Proportional width for encoding subplot
    trial_width = trial_data_size / total_size; % Proportional width for whole trial subplot
    
    % First subplot: Display window1_mean_PS_vectors (Item Encoding)
    subplot(2, 1, 1); % Create a 2x1 grid, first plot
    imagesc(ES_file.avg_window1_mean_PS_vectors');
    c=colorbar;
    c.Label.String = 'Z-scored Power';
    title('Item Encoding', 'Interpreter', 'none');
    xlabel('Encoding Windows');
    ylabel('Frequency');
    
    % Capture color limits from first plot, increase range and use for
    % both plots
    limits_to_use = clim;
    limits_to_use(1) = limits_to_use(1) - 1;
    limits_to_use(2) = limits_to_use(2) + 1;
    clim(limits_to_use)

    % Adjust subplot width based on data size ratio
    ax1 = gca; % Get current axis
    pos1 = get(ax1, 'Position'); % Get current position
    pos1(3) = encoding_width * 0.85; % Adjust width (0.85 to keep margins)
    set(ax1, 'Position', pos1); % Apply new position

    % Second subplot: Display window2_mean_PS_vectors (Whole Trial)
    subplot(2, 1, 2); % Create a 2x1 grid, second plot
    imagesc(ES_file.avg_window2_mean_PS_vectors');
    c=colorbar;
    clim(limits_to_use); % Apply same color limits as the first plot
    c.Label.String = 'Z-scored Power'; % Customize the label text as needed
    title('Whole Trial');
    xlabel('Trial Windows');
    ylabel('Frequency');
    
    % Plot the vertical lines
    index_adjust = -5; % subtracting 5 actually picks the window id that is centered at the desired time instead of the window id that begins at the desired time.
    xline(25+index_adjust, 'g', 'LineWidth', 2);    % baseline start
    xline(75+index_adjust, 'g', 'LineWidth', 2);    % baseline end
    xline(100+index_adjust, 'b', 'LineWidth', 2);   % fixation end
    xline(150+index_adjust, 'b', 'LineWidth', 2);   % enc1 end
    xline(200+index_adjust, 'b', 'LineWidth', 2);   % enc2 end
    xline(250+index_adjust, 'b', 'LineWidth', 2);   % enc3 end
    xline(650+index_adjust, 'b', 'LineWidth', 2);   % maintenance end
    % Get the upper limit of the y-axis for label placement
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

    % Adjust subplot width based on data size ratio
    ax2 = gca; % Get current axis
    pos2 = get(ax2, 'Position'); % Get current position
    pos2(3) = trial_width * 0.85; % Adjust width proportional to data size (0.85 to keep margins)
    set(ax2, 'Position', pos2); % Apply new position

    % Set overall figure title
    sgtitle(title_str, 'Interpreter', 'none'); % Set super title for the figure
    
    % Optionally, set colormap for both subplots
    % colormap('hot'); % Apply a colormap
end