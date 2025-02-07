function channels_to_process = subset_channels_by_brain_location(channel_IDs, brain_anatomies_to_process, patient_channel_brain_locations)
    %  brain_locations_to_process is a cell array of strings to match against
    %  patient_channel_brain_locations is a cell array of length original_N_channels
    %  channel_IDs are candidate IDs of the original data matrix (value corresponds to index in original matrix)
    if isempty(brain_anatomies_to_process)
        channels_to_process = channel_IDs;
    else
        % Initialize an empty logical array to store the indices to keep
        keep_indices = false(size(channel_IDs));
        
        % Loop through each channel in channels_to_process
        for i = 1:length(channel_IDs)
            channel_index = channel_IDs(i);
            
            % Get the brain location of the current channel (as string)
            channel_brain_location = patient_channel_brain_locations{channel_index};
            
            % Compare it with the brain locations to process
            for j = 1:length(brain_anatomies_to_process)
                % Check if the current brain location is partially in the brain location to process (case-insensitive)
                if contains(lower(channel_brain_location), lower(brain_anatomies_to_process{j}))
                    keep_indices(i) = true;  % Mark this index to keep
                    break;  % No need to check further if we already found a match
                end
            end
        end
        
        % Update channels_to_process to only include the ones that match
        channels_to_process = channel_IDs(keep_indices);
    
        % If no channels will be processed, throw an error
        if isempty(channels_to_process)
            warning('channels_to_process is empty. brain_anatomies_to_process must contain something relevant.');
        end
    end
end