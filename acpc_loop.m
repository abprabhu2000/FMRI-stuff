% Define the parent directory where all subject folders are located
parentDirectory = ['E:\flanker]flanker_task_test2'];  % Update this to your parent directory path

% Get a list of all subject directories within the parent directory
subjectDirs = dir(parentDirectory);

% Loop through all directories to find subjects and check their anat directories
for i = 1:length(subjectDirs)
    if subjectDirs(i).isdir && ~ismember(subjectDirs(i).name, {'.', '..'})  % Skip . and .. directories
        subjectName = subjectDirs(i).name;
        anatDir = fullfile(parentDirectory, subjectName, 'anat');  % Construct the path to the "anat" directory
        
        if exist(anatDir, 'dir')  % Check if the "anat" directory exists
            % Check if the anat folder already contains a "_T1w_reorint.mat" file
            % reorientedMatFile = fullfile(anatDir, sprintf('%s_T1w_reorint.mat', subjectName));
            
            % if exist(reorientedMatFile, 'file')
              %  fprintf('Skipping subject %s - ACPC correction already done.\n', subjectName);
               % continue;  % Skip this subject
            %end
            
            % Find all NIfTI files in the anat directory
            niftiFiles = dir(fullfile(anatDir, '*.nii'));
            
            if length(niftiFiles) == 1  % Check if there is exactly one NIfTI file
                niftiFilePath = fullfile(anatDir, niftiFiles(1).name);
                
                % Display the NIfTI file using SPM and wait until the window is closed
                fprintf('Loading NIfTI file for subject %s: %s\n', subjectName, niftiFilePath);
                spm_image('Display', niftiFilePath);
                
                % Wait for the user to close the SPM image window
                uiwait(msgbox(sprintf('Close the SPM window for subject %s to proceed to the next.', subjectName), 'Close SPM to Continue', 'modal'));
                
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

fprintf('Processing completed for all subjects.\n');
