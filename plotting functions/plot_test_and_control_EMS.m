function plot_test_and_control_EMS(chan_test_matrix, control_matrices, control_matrices_std, chan_id, chan_test_table, target_enc_ids, target_image_ids)
    % Compute the mean across dimensions after the first two
    chan_test_mean = mean(chan_test_matrix, 3);
    control_mean = mean(control_matrices, 3);
    control_std_mean = mean(control_matrices_std, 3);
    control_std_mean_after_averaging = std(control_matrices,0,3);
    

% for e=1:40
%     for m=1:390
%         p(e,m)=sum(chan_test_mean(e,m)>control_matrices(e,m,:))/1000;
%     end
% end
% 
% 
% figure; imagesc(p>0.9999999999)



    % Create a figure for the plots
    figure;
    clims=[0 0.5];
    % Plot the test matrix
    subplot(1, 4, 1);
    imagesc(chan_test_mean, clims);
    colorbar;
    title('Channel Test Mean');
    xlabel('Time');
    ylabel('Frequency');
    
    % Plot the mean of control matrices
    subplot(1, 4, 2);
    imagesc(control_mean, clims);
    colorbar;
    title('Mean of Control');
    xlabel('Time');
    ylabel('Frequency');
    
    % Plot the standard deviation of control matrices
    subplot(1, 4, 3);
    imagesc(control_std_mean, clims);
    colorbar;
    title('STD of Control (trial-level std)');
    xlabel('Time');
    ylabel('Frequency');

    % Plot the standard deviation of control matrices
    subplot(1, 4, 4);
    imagesc(control_std_mean_after_averaging, clims);
    colorbar;
    title('STD of Control (trial avg-level std)');
    xlabel('Time');
    ylabel('Frequency');
    
    % Add a super title with all the required details
    sgtitle(sprintf('Channel ID: %s, Anatomical Label: %s, Target Enc IDs: %s, Target Image IDs: %s', ...
        chan_id, chan_test_table.anatomical_label(1), mat2str(target_enc_ids), mat2str(target_image_ids)));


end