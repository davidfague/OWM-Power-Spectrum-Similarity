function channel_ids_to_keep = filter_channels_ids_to_use(channel_ids_to_use, original_labels)
    % Initialize the output array
    channel_ids_to_keep = [];
    
    % List of terms to exclude
    excluded_terms = ["pallidum", "putamen", "caudate", "ventricle", "white matter", "insula"];
    
    % Loop through the channel IDs
    for chan_idx = 1:length(channel_ids_to_use)
        chan_id = channel_ids_to_use(chan_idx);
        % Convert the label to lowercase for case-insensitive comparison
        label = lower(string(original_labels(chan_idx)));
        
        % Check if the label does not contain any of the excluded terms
        if all(~contains(label, excluded_terms))
            % Append the channel ID to the output array
            channel_ids_to_keep = [channel_ids_to_keep, chan_id];
        end
    end
end