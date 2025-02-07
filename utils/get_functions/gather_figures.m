% Specify the source folder (the one you want to duplicate)
sourceFolder = 'C:\Users\david\Downloads\Power Spectrum Similarity\output\';

% Specify the destination folder (the new location for only the .fig files)
destinationFolder = 'C:\Users\david\Downloads\Power Spectrum Similarity\whole trial P046 amyg enc3 correct\';

% Create the destination folder if it doesn't exist
if ~exist(destinationFolder, 'dir')
    mkdir(destinationFolder);
end

% Recursively get a list of all subfolders in the source folder
subfolders = genpath(sourceFolder);
subfolders = strsplit(subfolders, pathsep); % Split into individual folder paths

% Loop through each subfolder
for i = 1:length(subfolders)
    currentSubfolder = subfolders{i};
    
    if isempty(currentSubfolder)
        continue; % Skip empty folders
    end
    
    % Create the corresponding folder in the destination
    relativePath = strrep(currentSubfolder, sourceFolder, ''); % Relative path to preserve structure
    newSubfolder = fullfile(destinationFolder, relativePath);
    
    if ~exist(newSubfolder, 'dir')
        mkdir(newSubfolder);
    end
    
    % Get list of all .fig files in the current subfolder
    figFiles = dir(fullfile(currentSubfolder, '*.fig'));
    
    % Copy each .fig file to the new subfolder
    for j = 1:length(figFiles)
        sourceFile = fullfile(currentSubfolder, figFiles(j).name);
        destinationFile = fullfile(newSubfolder, figFiles(j).name);
        copyfile(sourceFile, destinationFile);
    end
end

disp('Folder structure copied with only .fig files.');
