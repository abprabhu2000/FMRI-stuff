
% Define the dictionary (equivalent to Python's `naming`)
naming = containers.Map({'ADHD', 'HC'}, {'5CSRTT-01', '5CSRTT-00'});

% Get a list of directories in the parent directory
dirs = dir(PD);

for i = 1:length(dirs)
    % Skip hidden/system files or non-directories
    if startsWith(dirs(i).name, '.') || startsWith(dirs(i).name, '_') || ~dirs(i).isdir
        continue;
    end

    new_filename = dirs(i).name;

    % Replace based on the naming dictionary
    keys = naming.keys;
    values = naming.values;

    for k = 1:length(keys)
        key = keys{k};
        value = values{k};
        if contains(new_filename, key)
            new_filename = strrep(new_filename, key, value);
        end
    end

    % If the new filename is different, rename the folder
    if ~strcmp(new_filename, dirs(i).name)
        new_file_path = fullfile(PD, new_filename);
        old_file_path = fullfile(PD, dirs(i).name);
        movefile(old_file_path, new_file_path);
        fprintf('%s renamed to %s\n', dirs(i).name, new_file_path);
    end
end
%% 
    % Define the dictionary (equivalent to Python's `naming`)
    naming = containers.Map({'ADHD', 'HC'}, {'5CSRTT-01', '5CSRTT-00'});
    
    % Get a list of subjects (subdirectories)
    subs = dir(PD);
    
    for i = 1:length(subs)
        % Skip hidden/system files or non-directories
        if startsWith(subs(i).name, '.') || startsWith(subs(i).name, '_') || ~subs(i).isdir
            continue;
        end
        
        subsp = fullfile(PD, subs(i).name);
        subdirs = dir(subsp);
        
        for j = 1:length(subdirs)
            % Skip hidden/system files or non-directories
            if startsWith(subdirs(j).name, '.') || startsWith(subdirs(j).name, '_') || ~subdirs(j).isdir
                continue;
            end
            
            slp = fullfile(subsp, subdirs(j).name);
            files = dir(slp);
            
            for k = 1:length(files)
                % Skip hidden/system files or non-files
                if startsWith(files(k).name, '.') || startsWith(files(k).name, '_') || files(k).isdir
                    continue;
                end
                
                fps = fullfile(slp, files(k).name);
                new_filename = files(k).name;
                
                % Replace based on the naming dictionary
                keys = naming.keys;
                values = naming.values;
                
                for m = 1:length(keys)
                    key = keys{m};
                    value = values{m};
                    if contains(new_filename, key)
                        new_filename = strrep(new_filename, key, value);
                    end
                end
                
                % If the new filename is different, rename the file
                if ~strcmp(new_filename, files(k).name)
                    new_file_path = fullfile(slp, new_filename);
                    movefile(fps, new_file_path);
                    fprintf('%s renamed to %s\n', fps, new_file_path);
                end
            end
        end
    end
%%

% moving 1st_Level
    % This function renames the '1st_Level' folder to the subject folder name
    % and moves the subject folder to the move_to directory.
    
    % Get list of all subjects in the parent directory (PD)
    subs = dir(PD);
    
    for i = 1:length(subs)
        subsp = fullfile(PD, subs(i).name);
        
        % Skip hidden/system files and non-directories
        if startsWith(subs(i).name, '.') || startsWith(subs(i).name, '_') || ~subs(i).isdir
            continue;
        end
        
        % List subfolders in the subject folder
        sublets = dir(subsp);
        
        for j = 1:length(sublets)
            slp = fullfile(subsp, sublets(j).name);
            
            % Skip hidden/system files and non-directories
            if startsWith(sublets(j).name, '.') || startsWith(sublets(j).name, '_') || ~sublets(j).isdir
                continue;
            end
            
            % If subfolder is '1st_Level', rename it to the subject name
            if strcmp(sublets(j).name, '1st_Level')
                new_folder_path = fullfile(subsp, subs(i).name);
                movefile(slp, new_folder_path);
                fprintf('%s/%s renamed to %s/%s\n', subs(i).name, sublets(j).name, subs(i).name, new_folder_path);
            end
            
            % If subfolder name matches subject name, move it to move_to
            if strcmp(sublets(j).name, subs(i).name)
                slp2 = fullfile(subsp, sublets(j).name);
                new_folder_path2 = fullfile(move_to, sublets(j).name);
                copyfile(slp2, new_folder_path2);
                fprintf('%s moved to %s\n', slp2, new_folder_path2);
            end
        end
    end

    %%

    % moving feedback 
    % This function renames the 'feedback' folder to the subject folder name
    % and moves the subject folder to the move_to directory.
    
    % Get list of all subjects in the parent directory (PD)
    subs = dir(PD);
    
    for i = 1:length(subs)
        subsp = fullfile(PD, subs(i).name);
        
        % Skip hidden/system files and non-directories
        if startsWith(subs(i).name, '.') || startsWith(subs(i).name, '_') || ~subs(i).isdir
            continue;
        end
        
        % List subfolders in the subject folder
        sublets = dir(subsp);
        
        for j = 1:length(sublets)
            slp = fullfile(subsp, sublets(j).name);
            
            % Skip hidden/system files and non-directories
            if startsWith(sublets(j).name, '.') || startsWith(sublets(j).name, '_') || ~sublets(j).isdir
                continue;
            end
            
            % If subfolder is 'feedback', rename it to the subject name
            if strcmp(sublets(j).name, 'feedback')
                new_folder_path = fullfile(subsp, subs(i).name);
                movefile(slp, new_folder_path);
                fprintf('%s/%s renamed to %s/%s\n', subs(i).name, sublets(j).name, subs(i).name, new_folder_path);
            end
            
            % If subfolder name matches subject name, move it to move_to
            if strcmp(sublets(j).name, subs(i).name)
                slp2 = fullfile(subsp, sublets(j).name);
                new_folder_path2 = fullfile(move_to, sublets(j).name);
                copyfile(slp2, new_folder_path2);
                fprintf('%s moved to %s\n', slp2, new_folder_path2);
            end
        end
    end

%%

feedback_dir = '';
niftis_dir = '';
    % This function goes through directories, finds 'rp_*.txt' files, 
    % and copies them to the 'func' directory within the corresponding folder in 'niftis_dir'.
    
    % Get a list of all subjects in the feedback directory
    subjects = dir(feedback_dir);
    
    % Loop through each subject
    for i = 1:length(subjects)
        subs = subjects(i).name;
        subsp = fullfile(feedback_dir, subs);
        
        % Skip hidden folders or non-directory items
        if startsWith(subs, '.') || startsWith(subs, '_') || ~isfolder(subsp)
            fprintf('%s is not a dir\n', subs);
            continue;
        end
        
        % Get the sub-directories and files in the current subject directory
        sublets = dir(subsp);
        
        % Loop through sub-directories or files
        for j = 1:length(sublets)
            sublName = sublets(j).name;
            slp = fullfile(subsp, sublName);
            
            % Skip hidden files or directories
            if startsWith(sublName, '.') || startsWith(sublName, '_') || isfolder(slp)
                continue;
            end
            
            % Check if the file starts with 'rp_' and ends with '.txt'
            if startsWith(sublName, 'rp_') && endsWith(sublName, '.txt')
                % Copy the file to the corresponding directory in 'niftis'
                newpath = fullfile(niftis_dir, subs, 'func', sublName);
                
                % Create the destination directory if it doesn't exist
                if ~exist(fullfile(niftis_dir, subs, 'func'), 'dir')
                    mkdir(fullfile(niftis_dir, subs, 'func'));
                end
                
                % Copy the file
                copyfile(slp, newpath);
                fprintf('Copied %s to %s\n', slp, newpath);
            end
        end
    end
    


