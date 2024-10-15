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
            % Path to the merged 4D NIfTI file
            oldFilePath = fullfile(funcDir, 'merged_4d_output.nii');
            
            % Check if the merged 4D NIfTI file exists
            if exist(oldFilePath, 'file')
                % Create the new filename based on the subject folder name
                newFileName = sprintf('%s_bold.nii', subjectName);  % For example: "HC-023_bold.nii"
                newFilePath = fullfile(funcDir, newFileName);
                
                % Rename the file
                try
                    movefile(oldFilePath, newFilePath);
                    fprintf('Renamed "%s" to "%s"\n', oldFilePath, newFilePath);
                catch ME
                    fprintf('Error renaming "%s": %s\n', oldFilePath, ME.message);
                end
            else
                fprintf('No "merged_4d_output.nii" file found in "%s"\n', funcDir);
            end
        else
            fprintf('No "func" directory found for subject: %s\n', subjectName);
        end
    end
end

fprintf('Renaming completed for all subjects.\n');
