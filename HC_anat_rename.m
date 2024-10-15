% Define the parent directory where all subject folders are located
parentDirectory = 'E:\data\subjects';  % Update this to your parent directory path

% Get a list of all subject directories within the parent directory
subjectDirs = dir(parentDirectory);

% Loop through all directories to find subjects containing "HC" in their name
for i = 1:length(subjectDirs)
    if subjectDirs(i).isdir && contains(subjectDirs(i).name, 'HC')  % Check if it is a directory and contains "HC"
        subjectName = subjectDirs(i).name;
        anatDir = fullfile(parentDirectory, subjectName, 'anat');  % Construct the path to the "anat" directory
        
        if exist(anatDir, 'dir')  % Check if the "anat" directory exists
            % Find all NIfTI files in the anat directory
            niftiFiles = dir(fullfile(anatDir, '*.nii'));
            
            if length(niftiFiles) == 1  % Check if there is exactly one NIfTI file
                oldFilePath = fullfile(anatDir, niftiFiles(1).name);
                
                % Create the new filename based on the subject folder name
                newFileName = sprintf('%s_T1w.nii', subjectName);  % For example: "HC-023_T1w.nii"
                newFilePath = fullfile(anatDir, newFileName);
                
                % Rename the file
                try
                    movefile(oldFilePath, newFilePath);
                    fprintf('Renamed "%s" to "%s"\n', oldFilePath, newFilePath);
                catch ME
                    fprintf('Error renaming "%s": %s\n', oldFilePath, ME.message);
                end
            elseif isempty(niftiFiles)
                fprintf('No NIfTI file found in "%s"\n', anatDir);
            else
                fprintf('Multiple NIfTI files found in "%s". Skipping...\n', anatDir);
            end
        else
            fprintf('No "anat" directory found for subject: %s\n', subjectName);
        end
    end
end

fprintf('Renaming completed for all subjects.\n');

