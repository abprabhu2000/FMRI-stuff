function renameFolders(parentDir)
    % Function to rename all folders named '1st_Level' to 'firstlevel' in the directory tree starting from parentDir.
    %
    % Parameters:
    % - parentDir (string): The parent directory to start the search from.
    %
    % Usage:
    % renameFolders('E:\your\parent\directory');

    % Get a list of all subdirectories in the parent directory
    dirList = dir(fullfile(parentDir, '**', '1st_Level'));  % Find all folders named '1st_Level'

    % Iterate through each found folder and rename it
    for i = 1:length(dirList)
        % Get the old folder name
        oldFolder = fullfile(dirList(i).folder, dirList(i).name);
        
        % Define the new folder name
        newFolder = fullfile(dirList(i).folder, 'firstlevel');

        % Rename the folder
        try
            movefile(oldFolder, newFolder);
            fprintf('Renamed: %s to %s\n', oldFolder, newFolder);
        catch ME
            fprintf('Failed to rename: %s\nError: %s\n', oldFolder, ME.message);
        end
    end
end
