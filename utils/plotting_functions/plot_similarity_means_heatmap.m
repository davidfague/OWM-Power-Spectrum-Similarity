function fig = plot_similarity_means_heatmap(WI,BI, plot_params)
        if plot_params.same_n
            BI = BI(:,:,randi(size(BI,3),size(WI,3), 1));
        end
        fig = figure;
        subplot(1,3,1)
        imagesc(mean(WI,3))
        colormap('jet')
        colorbar;
        clim([-0.1, 0.5])% clim([-0.15, 0.3])
        title(sprintf('Within-Images mean N=%s(trial pairs)', num2str(size(WI,3))))

        subplot(1,3,2)
        imagesc(mean(BI,3))
                colormap('jet')
        colorbar;
        clim([-0.1, 0.5])% clim([-0.15, 0.3])
        title(sprintf('Between-Images mean N=%s(trial pairs)', num2str(size(BI,3))))

        subplot(1,3,3)
        imagesc(mean(WI,3)-mean(BI,3))
        colormap('jet')
        colorbar;
        % clim([0, 0.3])% clim([-0.15, 0.3])
        title(sprintf('mean WI - mean BI'))


        title_str = sprintf("%s p%s chan%s image%s enc%s\n %s", ...
                plot_params.type, num2str(plot_params.patient_id), ... 
                num2str(plot_params.chan_id), num2str(plot_params.image_id), ...
                num2str(plot_params.enc_id), plot_params.anat);
        sgtitle(title_str);

        hold off;
end