function generate_dummy_nifti_files(outputFolder)
    % Function to generate 425 dummy NIfTI files with a specific naming convention.
    %
    % Parameters:
    %   outputFolder (string): Path to the folder where dummy NIfTI files will be saved.
    %
    % Example usage:
    %   generate_dummy_nifti_files('C:/data/nifti_files');
    
    % Check if output folder exists, if not, create it
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Define the number of dummy NIfTI files to generate
    numFiles = 425;
    
    % Initialize SPM defaults
    spm('defaults', 'FMRI');
    
    % Loop to generate dummy files
    for i = 1:numFiles
        % Create random 3D data (e.g., 64x64x36 volume)
        dummyData = rand(64, 64, 36);
        
        % Define the NIfTI file header
        niftiHeader = struct();
        niftiHeader.fname = fullfile(outputFolder, sprintf('%s-%05d-%05d-%06d-01.nii', ...
            generateRandomString(33), 12345, i, randi([100000, 999999])));
        niftiHeader.dim = size(dummyData);
        niftiHeader.dt = [16 0];  % Data type: 16 -> float32
        niftiHeader.mat = eye(4);  % Identity matrix for affine transformation
        niftiHeader.pinfo = [1; 0; 0];  % Scale and offset information
        
        % Write the dummy NIfTI file
        spm_write_vol(niftiHeader, dummyData);
    end
    
    fprintf('Generated %d dummy NIfTI files in folder: %s\n', numFiles, outputFolder);
end

function randStr = generateRandomString(strLength)
    % Function to generate a random alphanumeric string of specified length
    symbols = ['A':'Z' 'a':'z' '0':'9'];
    nums = randi(numel(symbols), [1 strLength]);
    randStr = symbols(nums);
end
