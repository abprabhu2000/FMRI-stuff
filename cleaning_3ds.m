% Define the parent directory where all subject folders are located
parentDirectory = 'E:\data\subjects';  % Update this to your parent directory path

% Get a list of all subject directories within the parent directory
subjectDirs = dir(parentDirectory);

% Loop through all directories to find subjects containing "HC" in their name
for i = 1:length(subjectDirs)
    if subjectDirs(i).isdir && contains(subjectDirs(i).name, 'HC')  % Check if it is a directory and contains "HC"
        subjectName = subjectDirs(i).name;
        funcDir = fullfile(parentDirectory, subjectName, 'func');  % Construct the path to the "func" directory

        if exist(funcDir, 'dir')  % Check if the "func" directory exists
            % Get a list of all files in the "func" directory
            funcFiles = dir(funcDir);

            % Loop through all files in the "func" directory
            for j = 1:length(funcFiles)
                fileName = funcFiles(j).name;
                filePath = fullfile(funcDir, fileName);

                % Skip directories and hidden files
                if funcFiles(j).isdir || startsWith(fileName, '.')
                    continue;
                end

                % Delete files except 'merged_4d_output.nii'
                if ~strcmp(fileName, 'merged_4d_output.nii')
                    try
                        delete(filePath);
                        fprintf('Deleted: %s\n', filePath);
                    catch ME
                        fprintf('Error deleting %s: %s\n', filePath, ME.message);
                    end
                end
            end
        else
            fprintf('No "func" directory found for subject: %s\n', subjectName);
        end
    end
end

fprintf('Cleanup completed for all subjects.\n');
