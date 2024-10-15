function merge_3d_to_4d_nii(inputFolder, outputFilename)
    % Function to merge 3D NIfTI files in a folder into a single 4D NIfTI file.
    % The merging is based on the order determined by the integer value of digits
    % at positions 41 to 45 in filenames.
    %
    % Parameters:
    %   inputFolder (string): Path to the folder containing 3D NIfTI files.
    %   outputFilename (string): Name of the output 4D NIfTI file.
    %
    % Example usage:
    %   merge_3d_to_4d_nii('C:/data/nifti_files', 'output_4d.nii');
    
    % Get list of all NIfTI files in the input folder
    niftiFiles = dir(fullfile(inputFolder, '*.nii'));
    
    % Extract the order number from characters 41 to 45 of each filename and convert to integer
    fileOrder = zeros(length(niftiFiles), 1);
    for i = 1:length(niftiFiles)
        filename = niftiFiles(i).name;
        fileOrder(i) = str2double(filename(41:45));  % Convert characters 41:45 to an integer
    end
    
    % Sort files based on the extracted integer order
    [~, idx] = sort(fileOrder);
    niftiFiles = niftiFiles(idx);
    
    % Initialize SPM defaults
    spm('defaults', 'FMRI');
    
    % Read all 3D NIfTI files and write them as a 4D file
    for i = 1:length(niftiFiles)
        % Load NIfTI file header and data
        niftiFilePath = fullfile(inputFolder, niftiFiles(i).name);
        niftiHeader = spm_vol(niftiFilePath);
        niftiData = spm_read_vols(niftiHeader);
        
        % Initialize 4D NIfTI on the first iteration
        if i == 1
            [x, y, z] = size(niftiData);
            mergedHeader = niftiHeader;
            mergedHeader.fname = outputFilename;
            mergedHeader.dim = [x, y, z];  % 3D dimension
            mergedHeader.dt = niftiHeader.dt;  % Use the data type of the first header
            mergedHeader.pinfo = [1; 0; 0];  % Standard scaling info for 4D file
            mergedHeader.n = [1, 1];  % Initialize volume index
            
            mergedHeader = spm_create_vol(mergedHeader);  % Create the first volume
        end
        
        % Update the volume index and write the 3D volume data into the 4D NIfTI file
        mergedHeader.n = [i, 1];
        spm_write_vol(mergedHeader, niftiData);
    end
    
    fprintf('Successfully merged 3D NIfTI files into 4D NIfTI file: %s\n', outputFilename);
end
