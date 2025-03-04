function common_channels = get_channels_from_both_channel_maps(WI_channel_map, BI_channel_map)
    % Get the keys for each channel map
    wi_keys = double(string(keys(WI_channel_map)));
    bi_keys = double(string(keys(BI_channel_map)));
    
    % Find the keys common to both maps
    common_channels = intersect(wi_keys, bi_keys);
end