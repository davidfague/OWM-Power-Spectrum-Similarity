    function patient_image_to_trialID_encodingID_encodingCorrectness = subset_by_encodingCorrectness_and_ID(patient_image_ids_to_keep, patient_image_to_trialID_encodingID_encodingCorrectness, target_correctness, target_encodingIDs)
        for image_idx = 1:length(patient_image_ids_to_keep)
            image_trialIDs = patient_image_to_trialID_encodingID_encodingCorrectness{image_idx, 1};
            image_encodingIDs = patient_image_to_trialID_encodingID_encodingCorrectness{image_idx, 2};
            image_enc_correctness = patient_image_to_trialID_encodingID_encodingCorrectness{image_idx, 3};
            [image_trialIDs, image_encodingIDs, image_enc_correctness] = get_image_trialIDs_by_encCorrectness(image_trialIDs, image_encodingIDs, image_enc_correctness, target_correctness, target_encodingIDs);
            patient_image_to_trialID_encodingID_encodingCorrectness{image_idx, 1} = image_trialIDs;
            patient_image_to_trialID_encodingID_encodingCorrectness{image_idx, 2} = image_encodingIDs;
            patient_image_to_trialID_encodingID_encodingCorrectness{image_idx, 3} = image_enc_correctness;
        end
    end