function fig = plot_WI_BI_significant_clusters(BI,WI, params, clusters_info)
        significant_clusters = [clusters_info.significant_positive_clusters, clusters_info.significant_negative_clusters];
        all_real_clusters = [clusters_info.real_positive_clusters, clusters_info.real_negative_clusters];
        % if params.same_n
        %     BI = BI(:,:,randi(size(BI,3),size(WI,3), 1));
        % end
        fig = figure('WindowState','maximized');
        subplot(1,3,1)
        imagesc(mean(WI,3))
        colormap('jet')
        colorbar;
        clim([0, 0.3])% clim([-0.15, 0.3])
        title(sprintf('Within-Images mean N=%s(trial pairs)', num2str(size(WI,3))))

        subplot(1,3,2)
        imagesc(mean(BI,3))
                colormap('jet')
        colorbar;
        clim([0, 0.3])% clim([-0.15, 0.3])
        title(sprintf('Between-Images mean N=%s(trial pairs)', num2str(size(BI,3))))
        
        subplot(1,3,3)
        is_sig_cluster = zeros(size(WI,[1 2]));
        if sum(significant_clusters) > 1
            is_sig_cluster(all_real_clusters(1).PixelIdxList{significant_clusters}) = 1;
        end
        imagesc(is_sig_cluster)
        title('Significant Clusters')

        title_str = sprintf("%s p%s chan%s image%s enc%s\n %s", ...
                params.type, num2str(params.patient_id), ... 
                num2str(params.chan_id), num2str(params.image_id), ...
                num2str(params.enc_id), params.anat);
        sgtitle(title_str);

        hold off;
end