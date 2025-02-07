function [image_ids_to_keep, image_to_trialID_encodingID_encodingCorrectness] = subset_trials_by_image_occurence(trial_IDs, images_to_process, enc_to_imageID, images, enc_correctness)
    % Subset enc_to_imageID and enc_correctness by trial_IDs
    [enc_to_imageID, enc_correctness] = subset_by_trial_IDs(trial_IDs, enc_to_imageID, enc_correctness);

    % Get a boolean array indicating which images should be kept with leniency
    keep_images = get_keep_images(images_to_process, images);
    if sum(keep_images) == 0 % no images found
        error('No patient images are among the selected images_to_process.')
    end
    % get array of integer image ids
    image_ids_to_keep = find(keep_images==1);

    % get the trials to subset and enc periods to match with this trial
    % subset.
    image_to_trialID_encodingID_encodingCorrectness = cell(length(image_ids_to_keep), 3);
    for idx = 1:length(image_ids_to_keep)
        image_id = image_ids_to_keep(idx);
        
        % Find the trial and encoding IDs corresponding to the image_id
        [image_trialIDs, image_encodingIDs] = find(enc_to_imageID == image_id);
        
        % Get the correctness for these trials and encodings
        % Use sub2ind to convert the row and column indices into linear indices
        enc_indices = sub2ind(size(enc_correctness), image_trialIDs, image_encodingIDs);
        % Extract the corresponding enc_correctness values using linear indices
        image_enc_correctness = enc_correctness(enc_indices);

        % Store the trial IDs, encoding periods, and correctness
        image_to_trialID_encodingID_encodingCorrectness{idx, 1} = image_trialIDs;      % Store trial IDs
        image_to_trialID_encodingID_encodingCorrectness{idx, 2} = image_encodingIDs;   % Store encoding periods
        image_to_trialID_encodingID_encodingCorrectness{idx, 3} = image_enc_correctness; % Store correctness
    end
    % clear idx image_id image_trialIDs image_encodingIDs image_enc_correctness

    % return image_ids_to_keep, image_to_trialID_encodingID_encodingCorrectness
end

function [enc_to_imageID_subset, enc_correctness_subset] = subset_by_trial_IDs(trial_IDs, enc_to_imageID, enc_correctness) % in case trial ID is already a subset
    enc_to_imageID_subset = enc_to_imageID(trial_IDs, :);
    enc_correctness_subset = enc_correctness(trial_IDs, :);
end

function keep_images = get_keep_images(images_to_process, images)
    if ~isempty(images_to_process)
        keep_images = false(size(images));
        images_lower = lower(images); % Convert to lowercase for case-insensitive comparison
        images_to_process_lower = lower(images_to_process);
    
        for i = 1:length(images_to_process_lower)
            % Use contains to allow lenient partial matching
            for j = 1:length(images_lower)
                if contains(images_lower{j}, images_to_process_lower{i}, 'IgnoreCase', true)
                    keep_images(j) = true;
                end
            end
        end
    else
        keep_images = true(size(images));
    end
end