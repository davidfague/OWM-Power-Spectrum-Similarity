function add_event_lines()
    % Plot vertical lines to indicate different phases
    index_adjust = -5; % Center adjustment
    xline(25 + index_adjust, 'g', 'LineWidth', 2);    % baseline start
    xline(75 + index_adjust, 'g', 'LineWidth', 2);    % baseline end
    xline(100 + index_adjust, 'b', 'LineWidth', 2);   % fixation end
    xline(150 + index_adjust, 'b', 'LineWidth', 2);   % enc1 end
    xline(200 + index_adjust, 'b', 'LineWidth', 2);   % enc2 end
    xline(250 + index_adjust, 'b', 'LineWidth', 2);   % enc3 end
    xline(650 + index_adjust, 'b', 'LineWidth', 2);   % maintenance end

    % Add text labels to vertical lines
    label_loc = ylim;
    label_loc = label_loc(2) + 5;
    add_text_label(25 + index_adjust, label_loc, 'baseline start', 'g');
    add_text_label(75 + index_adjust, label_loc, 'baseline end', 'g');
    add_text_label(100 + index_adjust, label_loc, 'fixation end', 'b');
    add_text_label(150 + index_adjust, label_loc, 'enc1 end', 'b');
    add_text_label(200 + index_adjust, label_loc, 'enc2 end', 'b');
    add_text_label(250 + index_adjust, label_loc, 'enc3 end', 'b');
    add_text_label(650 + index_adjust, label_loc, 'maint end', 'b');
end