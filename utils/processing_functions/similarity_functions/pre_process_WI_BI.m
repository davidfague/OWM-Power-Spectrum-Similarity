function [EMS_WI,EMS_BI] = pre_process_WI_BI(params, WI, BI) 
        % get rid of trials that were not all3 correct (only based on test
        % trials since control trials have already been filtered.)
        if isfield(WI, 'BT_ES')
            field_to_use = 'BT_ES'; % old files implementation
        else
            field_to_use = 'matrix'; % new channelmap implementation
        end
        if params.only_all3_correct
            rows_to_keep = sum(WI.chan_test_table.encoding_correctness,2) == 3;
            WI.(field_to_use) = WI.(field_to_use)(:, :, rows_to_keep, rows_to_keep);
            BI.(field_to_use) = BI.(field_to_use)(:,:,rows_to_keep,:);
        end

        if numel(BI.(field_to_use)) == 0 || numel(WI.(field_to_use)) == 0
            EMS_WI = WI.(field_to_use);
            EMS_BI = BI.(field_to_use);
            return
        end

        % Create logical indices to omit same-trial pairs for WI
        [~, logicalIdx2] = createLogicalIndex(size(WI.(field_to_use), 3));
    
        % Process EMS data
        EMS_WI = WI.(field_to_use)(:, :, logicalIdx2); % omits diagonal (same-trial pairs)
        EMS_BI = BI.(field_to_use)(:, :, :); % flatten (enc, maint) to (trial_pairs)
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
