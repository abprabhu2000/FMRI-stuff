function move_rp_files(parent_dir)
    % Function to move files with 'rp' in their name and '.txt' ending
    % from the 'func' directory to the '1st_Level' directory for each subject.
    %
    % Input:
    %   - parent_dir: Path to the parent directory containing subject folders.
    %
    % Example usage:
    %   move_rp_files('/data/parent_directory');

    % Get list of all subjects in the parent directory
    subjects = dir(parent_dir);
    subjects = subjects([subjects.isdir]);  % Keep only directories
    subjects = subjects(~ismember({subjects.name}, {'.', '..'}));  % Remove '.' and '..'

    % Loop through each subject directory
    for i = 1:length(subjects)
        subject_dir = fullfile(parent_dir, subjects(i).name);
        
        % Define 'func' and '1st_Level' directories for the subject
        func_dir = fullfile(subject_dir, 'func');
        level1_dir = fullfile(subject_dir, '1st_Level');

        % Check if 'func' directory exists
        if ~exist(func_dir, 'dir')
            warning('Directory %s does not exist. Skipping this subject.', func_dir);
            continue;
        end
        
        % Create '1st_Level' directory if it does not exist
        if ~exist(level1_dir, 'dir')
            mkdir(level1_dir);
        end

        % Find the file with 'rp' in its name and '.txt' ending
        rp_file = dir(fullfile(func_dir, 'rp*.txt'));
        
        if isempty(rp_file)
            warning('No rp*.txt file found in %s. Skipping this subject.', func_dir);
            continue;
        end
        
        % Move the file to the '1st_Level' directory
        try
            movefile(fullfile(func_dir, rp_file(1).name), fullfile(level1_dir, rp_file(1).name));
            fprintf('Moved file %s to %s\n', rp_file(1).name, level1_dir);
        catch ME
            warning('Could not move file %s: %s', rp_file(1).name, ME.message);
        end
    end
end
