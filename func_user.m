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
            % Specify the output filename for the merged 4D NIfTI file
            outputFilename = fullfile(funcDir, 'merged_4d_output.nii');

            % Print the current processing subject
            fprintf('Processing subject: %s\n', subjectName);

            % Try to merge the 3D NIfTI files into a 4D NIfTI file
            try
                merge_3d_to_4d_nii(funcDir, outputFilename);
                fprintf('Successfully merged files for subject: %s\n', subjectName);
            catch ME
                fprintf('Error processing subject %s: %s\n', subjectName, ME.message);
            end
        else
            fprintf('No "func" directory found for subject: %s\n', subjectName);
        end
    end
end

fprintf('Processing completed for all subjects.\n');

