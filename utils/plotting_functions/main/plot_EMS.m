function plot_EMS(WI, BI, params, folder_to_save_in)
    % plot EMS timecourse
    fig = plot_similarity_means_timecourse(WI,BI, params);
    saveas(fig, sprintf("%s/WI_BI_time_courses_same_n.fig", folder_to_save_in))
    
    % plot EMS heatmap
    for same_n = [true, false]
        params.same_n = same_n;
        fig = plot_similarity_means_heatmap(WI,BI, params);
        if params.same_n
            saveas(fig, sprintf("%s/WI_BI_heatmaps_same_n.fig", folder_to_save_in))
        else
            saveas(fig, sprintf("%s/WI_BI_heatmaps.fig", folder_to_save_in))
        end
    end
end