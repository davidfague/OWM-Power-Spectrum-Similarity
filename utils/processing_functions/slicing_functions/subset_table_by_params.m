function [rows_to_use] = subset_table_by_params(plot_params, label_table, without_image)
    % subset table by plot_params
    if ~without_image 
        % rows with the image and plot_params
        rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
        label_table.patient_ID == plot_params.patient_id & ... % patient
        label_table.encID_to_imageID(:,plot_params.enc_id)==plot_params.image_id & ... % image in first encoding
        sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
    else
        % rows without the image
        rows_to_use = label_table.channel_ID == plot_params.chan_id & ... % channel
        label_table.patient_ID == plot_params.patient_id & ... % patient
        sum(label_table.encID_to_imageID(:,:)~=plot_params.image_id,2)==3 & ... % image in first encoding
        sum(label_table.encoding_correctness(:,:), 2)==3; % all 3 correct
    end

end