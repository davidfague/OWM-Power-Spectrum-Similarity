function fig = generate_avg_EMS_vs_EFS_figure(EFS_file, EMS_file, folder_name)
    % can turn this into function and clean up
    fig = figure('WindowState','maximized', 'Visible', 'off');%, 'Visible', 'off');
    title_str = strrep(folder_name, '_', ' ');
    title_str = strcat('Average Within-Trial Encoding Similarity ', title_str);
    title(title_str)
    plot(EFS_file.average_EFS_vector, 'color','r')
    hold on
    plot(EMS_file.average_EMS_vector, 'color','b')
    yline(EFS_file.average_EFS, 'color','red')
    yline(EMS_file.average_EMS, 'color','blue')
    ylabel('Similarity (rho)')
    xlabel('time window index')
    legend('enc-fixation', 'enc-maint', 'average EFS','average EMS')
end