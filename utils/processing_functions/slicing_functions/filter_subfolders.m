function [subfolders_to_analyze] = filter_subfolders(subFolders, patients_to_analyze, brain_locations_to_analyze, channels_to_analyze, gamma_chans_to_analyze, images_to_analyze, trials_to_analyze, encID_to_analyze, encCorrectness_to_analyze)
    % Initialize an array to keep track of which folders to keep
    keep_folders = false(length(subFolders), 1); 
    % Loop over each folder
    for k = 1:length(subFolders)
        name = subFolders(k).name;
        parts = strsplit(name, '_');  % Split the folder name into parts
        
        % Extract parts for comparison
        patientPart = str2double(regexprep(parts{1}, 'Patient', ''));
        brainLocationPart = parts{2};
        channelPart = regexprep(parts{3}, 'Chan', '');
        gammaChanPart = logical(str2double(regexprep(parts{4}, 'ChanisGam', '')));
        imagePart = regexprep(parts{5}, 'Image', '');
        trialPart = regexprep(parts{6}, 'Trial', '');
        encIdPart = str2double(regexprep(parts{7}, 'EncId', ''));
        correctnessPart = str2double(regexprep(parts{8}, 'Correct', ''));
        
        % Check if folder meets all criteria
        patientMatch = isempty(patients_to_analyze) || ismember(patientPart, patients_to_analyze);
        brainLocationMatch = isempty(brain_locations_to_analyze) || any(contains(lower(brainLocationPart), lower(brain_locations_to_analyze)));
        channelMatch = isempty(channels_to_analyze) || any(contains(channelPart, channels_to_analyze));
        gammaChanMatch = isempty(gamma_chans_to_analyze) || ismember(gammaChanPart, gamma_chans_to_analyze);
        imageMatch = isempty(images_to_analyze) || any(contains(imagePart, images_to_analyze));
        trialMatch = isempty(trials_to_analyze) || ismember(trialPart, trials_to_analyze);
        encIdMatch = isempty(encID_to_analyze) || ismember(encIdPart, encID_to_analyze);
        correctnessMatch = isempty(encCorrectness_to_analyze) || ismember(correctnessPart, encCorrectness_to_analyze);
        
        % If all criteria match, keep this folder
        if patientMatch && brainLocationMatch && channelMatch && gammaChanMatch && imageMatch && trialMatch && encIdMatch && correctnessMatch
            keep_folders(k) = true;
        end
    end
    
    % Filter the list of folders
    subfolders_to_analyze = subFolders(keep_folders);