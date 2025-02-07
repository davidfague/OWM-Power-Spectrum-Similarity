function channels_to_process = subset_channels_by_gamma_modulation(original_channel_IDs, use_gamma_mod_chans, is_gamma_channels)
    % original_channel_IDs: candidate IDs of the original data matrix (value corresponds to index in original matrix)
    % use_gamma_mod_chans: a logical array indicating whether to include gamma-modulated (true) and/or non-gamma-modulated (false) channels.
    % is_gamma_channels: a logical array indicating which channels are gamma modulated
    
    % Initialize channels_to_process as an empty array
    channels_to_process = [];
    
    % Check if original_channel_IDs is already a subset of channels
    % (subsets is_gamma_channels by channels_IDs in case channel IDs are a
    % subset of the original channels)
    candidate_is_gamma_channels = is_gamma_channels(original_channel_IDs);
    
    % If use_gamma_mod_chans contains true, include gamma-modulated channels
    if any(use_gamma_mod_chans == true)
        channels_to_process = [channels_to_process, original_channel_IDs(candidate_is_gamma_channels)];
    end
    
    % If use_gamma_mod_chans contains false, include non-gamma-modulated channels
    if any(use_gamma_mod_chans == false)
        channels_to_process = [channels_to_process, original_channel_IDs(~candidate_is_gamma_channels)];
    end
    
    % If no channels are selected, throw an error
    if isempty(channels_to_process)
        error('use_gamma_mod_chans must contain true, false, or both.');
    end
end
