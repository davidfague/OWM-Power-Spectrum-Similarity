function [EMS_WI,EMS_BI] = get_and_pre_process_WI_BI(patient_id, comp_options, enc_id, image_id, chan_id, params, type, only_all3_correct)
        if nargin < 8
            only_all3_correct = false;
        end

        [BI, WI] = load_BT_data(patient_id, comp_options, enc_id, image_id, chan_id, params, type); % doesn't matter that it looks like type='EMS' is returned.

        % get rid of trials that were not all3 correct (only based on test
        % trials since control trials have already been filtered.)
        if only_all3_correct
            rows_to_keep = sum(WI.chan_test_table.encoding_correctness,2) == 3;
            WI.BT_ES = WI.BT_ES(:, :, rows_to_keep, rows_to_keep);
            BI.BT_ES = BI.BT_ES(:,:,rows_to_keep,:);
        end

        if numel(BI.BT_ES) == 0 || numel(WI.BT_ES) == 0
            EMS_WI = WI.BT_ES;
            EMS_BI = BI.BT_ES;
            return
        end

        % Create logical indices to omit same-trial pairs for WI
        [logicalIdx1, logicalIdx2] = createLogicalIndex(size(WI.BT_ES, 3));
    
        % Process EMS data
        EMS_WI = WI.BT_ES(:, :, logicalIdx2); % omits diagonal (same-trial pairs)
        EMS_BI = BI.BT_ES(:, :, :);
        EMS_WI = EMS_WI(:, :, any(~isnan(EMS_WI), [1 2])); % Remove NaNs from WI
        EMS_BI = EMS_BI(:, :, any(~isnan(EMS_BI), [1 2])); % Remove NaNs from BI
end

function [logicalIdx1, logicalIdx2] = createLogicalIndex(sizeDim) % omits diagonal for WI pairs
    % createLogicalIndex creates two logical arrays for indexing a matrix,
    % excluding diagonal elements (e.g., 1,1; 2,2; ...; n,n).
    %
    % INPUT:
    %   sizeDim - Size of the dimension (e.g., 12 for a matrix (:,:,12,12))
    %
    % OUTPUT:
    %   logicalIdx1 - Logical index for the first dimension
    %   logicalIdx2 - Logical index for the second dimension

    % Validate input
    if ~isscalar(sizeDim) || sizeDim <= 0 || sizeDim ~= round(sizeDim)
        error('Input sizeDim must be a positive integer. %s', num2str(sizeDim));
    end

    % Initialize logical arrays
    logicalIdx1 = true(sizeDim, sizeDim);
    logicalIdx2 = true(sizeDim, sizeDim);

    % Create diagonal indices to exclude
    diagIdx = 1:sizeDim;

    % Set diagonal elements to false
    logicalIdx1(sub2ind([sizeDim, sizeDim], diagIdx, diagIdx)) = false;
    logicalIdx2(sub2ind([sizeDim, sizeDim], diagIdx, diagIdx)) = false;
end
