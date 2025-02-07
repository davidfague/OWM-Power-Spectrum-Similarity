function plot_test_and_control_EMS(chan_test_matrix, control_matrices, control_matrices_std, chan_id, chan_test_table, target_enc_ids, target_image_ids)
save=false;
% control_matrices 41, 390, 1, 1000
% chan_test_matrix 18, 41, 640
chan_test_matrix = chan_test_matrix(:, :, 251:640); % chan_test_matrix 18, 41, 390
chan_test_mean = squeeze(mean(chan_test_matrix, 1));

% % Compute the mean across dimensions after the first two
%     chan_test_mean = mean(chan_test_matrix, 3);
%     control_mean = mean(control_matrices, 3);
%     control_std_mean = mean(control_matrices_std, 3);
%     control_std_mean_after_averaging = std(control_matrices,0,3);
    

p = nan(size(chan_test_mean));
for e=1:41
    for m=1:390
        p(e,m)=sum(chan_test_mean(e,m)>control_matrices(e,m,:))/1000;
    end
end
% 
% 
% figure; imagesc(p>0.9999999999)


    % Create a figure for the plots
    figure;
    clims=[0 0.5];

    % Plot the test avg matrix
    subplot(3, 3, 1);
    imagesc(squeeze(chan_test_matrix(1,:,:)), clims);
    colorbar;
    title('Test - single trial');
    xlabel('Time');
    ylabel('Enc Time');

    % Plot the test avg matrix
    subplot(3, 3, 2);
    imagesc(chan_test_mean, clims);
    colorbar;
    title('Test - avg trial');
    xlabel('Time');
    ylabel('Enc Time');

    % Plot the test std
    subplot(3, 3, 3);
    imagesc(squeeze(std(chan_test_matrix, 0, 1)), clims);
    colorbar;
    title('Test - std across trial');
    xlabel('Time');
    ylabel('Enc Time');

    % Plot the single ctrl matrix
    subplot(3, 3, 4);
    imagesc(control_matrices(:,:,:,1), clims);
    colorbar;
    title('Control - single sample avg (single trialpair n.a.)');
    xlabel('Time');
    ylabel('Enc Time');

    % Plot the ctrl avg matrix
    subplot(3, 3, 5);
    imagesc(squeeze(mean(control_matrices, 4)), clims);
    colorbar;
    title('Control - sample avg (avg of avg across trial pairs');
    xlabel('Time');
    ylabel('Enc Time');

    % Plot the ctrl avg matrix
    subplot(3, 3, 6);
    imagesc(squeeze(mean(control_matrices_std, 3)), clims);
    colorbar;
    title('Control - std (avg std across trialpairs)');
    xlabel('Time');
    ylabel('Enc Time');

    % p values
    subplot(3, 3, 7);
    imagesc(p, [0 1]);
    colorbar;
    title('raw percentile of test mean on control sample means');
    xlabel('Time');
    ylabel('Enc Time');

    subplot(3, 3, 8);
    imagesc(p> .99, [0 1]);
    colorbar;
    title('percentile > .99');
    xlabel('Time');
    ylabel('Enc Time');

    % Plot the ctrl std across avg of trial pairs
    subplot(3, 3, 9);
    imagesc(squeeze(std(control_matrices, 0, 4)), clims);
    colorbar;
    title('Control - std (std across avg of trial pairs)');
    xlabel('Time');
    ylabel('Enc Time');

    %% old
    % % Create a figure for the plots
    % figure;
    % clims=[0 0.5];
    % % Plot the test matrix
    % subplot(1, 4, 1);
    % imagesc(chan_test_mean, clims);
    % colorbar;
    % title('Channel Test Mean');
    % xlabel('Time');
    % ylabel('Frequency');
    % 
    % % Plot the mean of control matrices
    % subplot(1, 4, 2);
    % imagesc(control_mean, clims);
    % colorbar;
    % title('Mean of Control');
    % xlabel('Time');
    % ylabel('Frequency');
    % 
    % % Plot the standard deviation of control matrices
    % subplot(1, 4, 3);
    % imagesc(control_std_mean, clims);
    % colorbar;
    % title('STD of Control (trial-level std)');
    % xlabel('Time');
    % ylabel('Frequency');
    % 
    % % Plot the standard deviation of control matrices
    % subplot(1, 4, 4);
    % imagesc(control_std_mean_after_averaging, clims);
    % colorbar;
    % title('STD of Control (trial avg-level std)');
    % xlabel('Time');
    % ylabel('Frequency');
    % 
    % % Add a super title with all the required details
    sgtitle(sprintf('Channel ID: %d, Anatomical Label: %s, Target Enc IDs: %s, Target Image IDs: %s', ...
        chan_id, chan_test_table.anatomical_label(1), mat2str(target_enc_ids), mat2str(target_image_ids)));

    if save
        saveas(figHandle, fullfile(saveDir, [figName, '.png']));
        closefig()
    end
end