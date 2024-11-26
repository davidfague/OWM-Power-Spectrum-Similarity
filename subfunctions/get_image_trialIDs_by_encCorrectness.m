function [selected_trialIDs, selected_image_encodingIDs, selected_image_enc_correctness] = get_image_trialIDs_by_encCorrectness(image_trialIDs, image_encodingIDs, image_enc_correctness, target_correctness, target_encodingIDs)
    % Find logical indices where the correctness matches and the encoding ID is in the target list
    valid_idx = (image_enc_correctness == target_correctness) & ismember(image_encodingIDs, target_encodingIDs);
    
    % Apply the logical index to select the corresponding values
    selected_trialIDs = image_trialIDs(valid_idx);
    selected_image_encodingIDs = image_encodingIDs(valid_idx);
    selected_image_enc_correctness = image_enc_correctness(valid_idx);
end
