% Specify your new base folder where all subfolders will be created
new_base_folder = 'C:\Users\david\Downloads\Power Spectrum Similarity\whole trial P046 amyg enc3 correct figures\';

% Get all open figures
figs = findall(0, 'Type', 'figure');

% Loop through all figures
for i = 1:length(figs)
    % Get the figure handle
    fig = figs(i);

    % Try to get the figure's name or number for generating a file name
    fig_name = get(fig, 'Name');  % Get the figure's 'Name' property (if set)
    
    if isempty(fig_name)
        % If the figure name is not set, use the figure number instead
        fig_name = ['Figure_' num2str(fig.Number)];
    end

    % Create a new subfolder based on the figure name
    new_subfolder = fullfile(new_base_folder, fig_name);

    % Create the new subfolder if it does not exist
    if ~exist(new_subfolder, 'dir')
        mkdir(new_subfolder);
    end

    % Construct the new filename for the figure in the new folder
    new_filename = fullfile(new_subfolder, [fig_name '.fig']);  % Save as .fig (MATLAB format)

    % Save the figure in the new subfolder
    % saveas(fig, new_filename);

    % Optionally save as other formats, such as PNG or JPEG
    saveas(fig, fullfile(new_subfolder, [fig_name '.png']));  % Save as PNG
end