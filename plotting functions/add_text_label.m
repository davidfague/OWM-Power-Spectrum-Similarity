function add_text_label(x, y, text_str, color)
    % Helper function to add text labels at specific x locations
    text(x, y, text_str, 'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', color);
end