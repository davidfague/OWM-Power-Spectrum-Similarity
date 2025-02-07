function fig = plot_PSV(whole_trial_PSV_data, plot_params, title_str, enforced_clim)
    if nargin< 4
        enforced_clim = [];
    end
    if ~isempty(plot_params.frequencies_to_use)
        whole_trial_PSV_data = whole_trial_PSV_data(:,plot_params.frequencies_to_use);
    end

    fig = figure('WindowState','maximized');%, 'Visible', 'off');%,'Visible', 'off');
    
    % title_str = sprintf("p%s chan%s %s image%s", num2str(plot_params.patient_id), ...
    %     num2str(plot_params.chan_id), plot_params.anat, num2str(plot_params.image_id));
    title_str = sprintf("Power Spectrum Vectors (PSVs)\n %s \n %s", title_str, plot_params.anat);

    % Calculate the relative widths based on the data sizes
    encoding_data_size = 40; % Encoding x-axis size (40)
    trial_data_size = size(whole_trial_PSV_data, 1); % Whole trial x-axis size (891)
    total_size = encoding_data_size + trial_data_size;
 
    encoding_width = encoding_data_size / total_size; % Proportional width for encoding subplot
    trial_width = trial_data_size / total_size; % Proportional width for whole trial subplot
    
    % First subplot: Display window1_mean_PS_vectors (Item Encoding)
    subplot(2, 1, 1); % Create a 2x1 grid, first plot
    imagesc(whole_trial_PSV_data(plot_params.enc_window_ids,:)');
    c=colorbar;
    c.Label.String = 'Z-scored Power';
    title('Item Encoding', 'Interpreter', 'none');
    xlabel('Encoding Windows');
    ylabel('Frequency');
    
    % Capture color limits from first plot, increase range and use for
    % both plots
    % limits_to_use = clim;
    % limits_to_use(1) = limits_to_use(1) - 1;
    % limits_to_use(2) = limits_to_use(2) + 1;
    % clim(limits_to_use)
    % Calculate the mean and standard deviation of the data to set clim
    % Flatten the data to a vector
    if isempty(enforced_clim)
        flattened_data = whole_trial_PSV_data(plot_params.enc_window_ids, :);
        % Calculate the mean and standard deviation of the flattened data
        data_mean = mean(flattened_data(:));
        data_std = std(flattened_data(:));
        % Set the color limits to be 2 standard deviations around the mean
        limits_to_use = [data_mean - 2*data_std, data_mean + 2*data_std];
    else
        limits_to_use = enforced_clim;
    end
    clim(limits_to_use)

    % Adjust subplot width based on data size ratio
    ax1 = gca; % Get current axis
    pos1 = get(ax1, 'Position'); % Get current position
    pos1(3) = encoding_width * 0.85; % Adjust width (0.85 to keep margins)
    set(ax1, 'Position', pos1); % Apply new position

    % Second subplot: Display window2_mean_PS_vectors (Whole Trial)
    subplot(2, 1, 2); % Create a 2x1 grid, second plot
    imagesc(whole_trial_PSV_data');
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