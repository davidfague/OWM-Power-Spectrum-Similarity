function fig = plot_similarity_means(WI,BI, plot_params)
        fig = figure;
        subplot(2,1,1); title('mean'); hold on; % means
        legend("AutoUpdate","on")
        % subplot(1,2,1)
        % imagesc(mean(WI,3))
        % clim([-0.15, 0.3])
        plot(mean(WI,[1 3]), 'DisplayName', 'Within-Images')
        % ylim([0 0.2])
        yline(mean(WI, [1:3]),'Label', 'Within-Images mean', 'DisplayName', 'Within-Images mean')
        % title(sprintf('WI mean %s', num2str(size(WI,3))))

        % subplot(1,2,2)
        % imagesc(mean(BI,3))
        % clim([-0.15, 0.3])
        plot(mean(BI,[1 3]), 'DisplayName', 'Between-Images')
        % ylim([0 0.2])
        yline(mean(BI, [1:3]), 'Label', 'Between-Images mean', 'DisplayName', 'Between-Images mean')
        % title(sprintf('BI mean %s', num2str(size(BI,3))))

        subplot(2,1,2); hold on;% stds
        legend("AutoUpdate","on")
        plot(std(WI,0,[1 3]), 'DisplayName', 'Within-Images'); hold on;
        % ylim([0 0.2])
        yline(std(WI, 0, [1:3]),'Label', 'Within-Images mean', 'DisplayName', 'Within-Images all obs')

        plot(std(BI,0,[1 3]), 'DisplayName', 'Between-Images')
        % ylim([0 0.2])
        yline(std(BI, 0, [1:3]), 'Label', 'Between-Images mean', 'DisplayName', 'Between-Images all obs')
        title('standard deviaton')

        title_str = sprintf("%s p%s chan%s image%s enc%s\n %s", ...
                plot_params.type, num2str(plot_params.patient_id), ... 
                num2str(plot_params.chan_id), num2str(plot_params.image_id), ...
                num2str(plot_params.enc_id), plot_params.anat);
        sgtitle(title_str);

        hold off;
end