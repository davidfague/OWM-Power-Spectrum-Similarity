function fig = plot_similarity_means_timecourse(WI,BI, plot_params)
        % Define colors (these are MATLAB's default "Lines" colors #1 and #2)
        WI_color = [0.0000, 0.4470, 0.7410];  % blue
        BI_color = [0.8500, 0.3250, 0.0980];  % reddish-orange

        mean_WI = squeeze(mean(WI, [1 3]));
        mean_BI = squeeze(mean(BI, [1 3]));
        WI_err = squeeze(std(mean(WI, 1), 0, [1 3]) / sqrt(size(WI, 3)));     % Standard deviation (0 flag for normalization by N-1) %  / sqrt(n) to convert to standard error of the mean
        BI_err = squeeze(std(mean(BI, 1), 0, [1 3]) / sqrt(size(WI, 3)));

        if ~strcmp(plot_params.type, 'EMS')
            error('%s not implemented in plot_similarity_means_timecourse', plot_params.type);
        else
            times = 1:length(mean_WI);
            times = (times + 260) * 10; % move to maintenance and convert to ms from 10 ms
        end

        % Calculate the upper and lower bounds of the error region
        mean_WI_upper = mean_WI + WI_err;
        mean_WI_lower = mean_WI - WI_err;
        mean_BI_upper = mean_BI + BI_err;
        mean_BI_lower = mean_BI - BI_err;

        fig = figure;
        subplot(4,1,1); title('mean'); hold on; % means
        legend("AutoUpdate","on")
        % subplot(1,2,1)
        % imagesc(mean(WI,3))
        % clim([-0.15, 0.3])
        % ylim([0 0.2])
        fill([times, fliplr(times)], [mean_WI_upper, fliplr(mean_WI_lower)], ...
            WI_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, ... % Make the fill semi-transparent
            'DisplayName', 'WI SEM');
        plot(times, mean_WI, 'DisplayName', 'Within-Images')

        yline(mean(WI, [1:3]),'Label', 'Within-Images overall mean', 'DisplayName', 'Within-Images overall mean')
        % title(sprintf('WI mean %s', num2str(size(WI,3))))

        % subplot(1,2,2)
        % imagesc(mean(BI,3))
        % clim([-0.15, 0.3])
        % --- Between-Images shaded area & line ---
        fill([times, fliplr(times)], ...
             [mean_BI_upper, fliplr(mean_BI_lower)], ...
             BI_color, ...
             'FaceAlpha', 0.2, ...
             'EdgeColor', 'none', ...
             'DisplayName', 'BI SEM');
        plot(times, mean_BI, 'DisplayName', 'Between-Images')
        % ylim([0 0.2])
        yline(mean(BI, [1:3]), 'Label', 'Between-Images overall mean', 'DisplayName', 'Between-Images overall mean')
        % title(sprintf('BI mean %s', num2str(size(BI,3))))

        subplot(4,1,2); % individual trials
        plot(times, squeeze(mean(WI,1)))
        title('WI trials')
        subplot(4,1,3); % individual trials
        plot(times, squeeze(mean(BI,1)))
        title('BI trials')

        subplot(4,1,4); hold on;% stds
        legend("AutoUpdate","on")
        plot(times, std(WI,0,[1 3]), 'DisplayName', 'Within-Images'); hold on;
        % ylim([0 0.2])
        yline(std(WI, 0, [1:3]),'Label', 'Within-Images mean', 'DisplayName', 'Within-Images all obs')

        plot(times,std(BI,0,[1 3]), 'DisplayName', 'Between-Images')
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