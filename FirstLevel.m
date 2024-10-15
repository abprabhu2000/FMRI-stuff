classdef FirstLevel
    % FirstLevel Class
    %
    % This class is designed to facilitate fMRI data analysis using SPM (Statistical
    % Parametric Mapping). It provides methods to prepare directories and files,
    % sort and move fMRI files, run SPM batch processing, and clean up the output
    % directories. The class handles all necessary configurations for an fMRI
    % analysis, including setting up model specifications, estimating the GLM, and
    % generating contrasts.
    %
    % **Order of Class Function Usage:**
    % 1. **filesorter:** Creates necessary directories within each subject's folder.
    % 2. **filemover:** Moves required fMRI files into the 'firstlevel' directory.
    % 3. **firstrun:** Runs the SPM batch processing for each subject's fMRI data.
    % 4. **FirstLevelCleaner:** Cleans up the output files by organizing them into
    %    respective folders based on file types.
    %
    % **Properties:**
    % - `PD` (Parent Directory): The root directory where the subject data is located.
    % - `RT` (Repetition Time): The duration of stimuli (in seconds).
    % - `Units` (Units of Time Measurement): Units of time measurement for stimuli duration.

    properties
        PD  % Parent directory
        RT  % Repetition Time
        Units  % Units of time measurement for duration of stimuli
    end

    methods
        % Constructor
        function obj = FirstLevel(parent, RT, units)
            % Constructor: FirstLevel
            %
            % Initializes an instance of the FirstLevel class with the specified parent
            % directory, repetition time, and time units.
            %
            % **Usage:**
            % ```
            % obj = FirstLevel(parent, RT, units)
            % ```
            %
            % **Parameters:**
            % - `parent` (char): Path to the parent directory containing subject folders.
            % - `RT` (double): Repetition time for the fMRI data (in seconds).
            % - `units` (char): Time units for stimuli durations (e.g., 'secs').
            %
            % **Returns:**
            % - `obj`: An instance of the FirstLevel class.

            obj.PD = parent;  % Parent directory
            obj.RT = RT;  % Duration of stimuli (Repetition Time)
            obj.Units = units;  % Units of time measurement
        end

        % Method to sort files and create necessary directories
        function filesorter(obj)
            % Method: filesorter
            %
            % Creates the necessary directories ('firstlevel', 'Res', 'Betas', 'Preproc',
            % 'Mats', 'Con', 'SPM_T') within each subject's folder for storing analysis
            % outputs and intermediary files.
            %
            % **Usage:**
            % ```
            % obj.filesorter()
            % ```
            %
            % **Parameters:**
            % - None
            %
            % **Returns:**
            % - None
            pd = obj.PD;
            subjects = dir(pd);

            for i = 1:length(subjects)
                subName = subjects(i).name;
                subf = fullfile(pd, subName);

                if ~isfolder(subf) || startsWith(subName, '.')
                    continue;
                end

                % Create 'firstlevel' directory and its subdirectories
                flp = fullfile(subf, 'firstlevel');
                if ~exist(flp, 'dir')
                    mkdir(flp);
                end

                % Create 'Res', 'Betas', 'Preproc', 'Mats' directories
                mkdir(fullfile(flp, 'Res'));
                mkdir(fullfile(flp, 'Betas'));
                mkdir(fullfile(flp, 'Preproc'));
                mkdir(fullfile(flp, 'Mats'));
                mkdir(fullfile(flp,'Con'));
                mkdir(fullfile(flp,'SPM_T'));

            end
        end

        % Method to move files to 'firstlevel' directory
        function filemover(obj,unique_term)
            % Method: filemover
            %
            % Moves functional fMRI files from the 'func' folder to the 'firstlevel' directory
            % within each subject's folder, based on a unique term that identifies the files.
            %
            % **Usage:**
            % ```
            % obj.filemover(unique_term)
            % ```
            %
            % **Parameters:**
            % - `unique_term` (string): A unique term to identify the fMRI files to be moved
            %   (e.g., 'swarADHD'). Only files containing this term and having a `.nii`
            %   extension will be considered for moving.
            %
            % **Behavior:**
            % - For each subject folder:
            %   1. Checks if the 'firstlevel' directory exists.
            %   2. Iterates through the files in the 'firstlevel' directory and skips files
            %      that already contain the `unique_term` and have the `.nii` extension.
            %   3. If such files are not found in 'firstlevel', searches in the 'func'
            %      directory within the subject folder.
            %   4. If matching files are found in 'func', they are moved to 'firstlevel'.
            %   5. If a file already exists in the destination ('firstlevel'), it is not
            %      copied again.
            %
            % **Returns:**
            % - None
            %
            % **Examples:**
            % ```
            % % Moves fMRI files containing the term 'swarADHD' from 'func' to 'firstlevel'
            % obj.filemover('swarADHD')
            % ```

            pd = obj.PD;
            subjects = dir(pd);

            for i = 1:length(subjects)
                subName = subjects(i).name;
                subf = fullfile(pd, subName);

                if ~isfolder(subf) || startsWith(subName, '.')
                    continue;
                end

                flp = fullfile(subf, 'firstlevel');
                if isfolder(flp)
                    files = dir(flp);
                    for j = 1:length(files)
                        fileName = files(j).name;
                        if contains(fileName, unique_term) && endsWith(fileName, '.nii') && isstring(unique_term)
                            continue;
                        else
                            req_file_folder = fullfile(subf, 'func');
                            if isfolder(req_file_folder)
                                funcs = dir(req_file_folder);
                                for k = 1:length(funcs)
                                    funcName = funcs(k).name;
                                    if contains(funcName, unique_term) && endsWith(funcName, '.nii') && isstring(unique_term)
                                        funcs_path_old = fullfile(req_file_folder, funcName);
                                        funcs_path_new = fullfile(flp, funcName);
                                        if ~exist(funcs_path_new, 'file')
                                            copyfile(funcs_path_old, funcs_path_new);
                                            fprintf('Moved %s to %s\n', funcs_path_old, funcs_path_new);
                                        else
                                            fprintf('File %s already exists.\n', funcs_path_new);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        % Method to run the SPM batch processing directly
        function firstRunAll(obj, NumConditions, cond_name_list, pmods,uqIdentifier)
            % Method: firstrun
            %
            % Runs SPM batch processing for each subject's fMRI data, including model
            % specification, GLM estimation, and contrast generation. The function creates
            % and runs a batch using the provided number of conditions, condition names,
            % and parametric modulation (pmod) details.
            %
            % **Usage:**
            % ```
            % obj.firstrun(NumConditions, cond_name_list, pmods)
            % ```
            %
            % **Parameters:**
            % - `NumConditions` (int): Number of conditions in the fMRI experiment.
            % - `cond_name_list` (cell array of char): List of condition names.
            % - `pmods` (cell array of char): List of conditions that have parametric modulation.
            %
            % **Returns:**
            % - None
            pd = obj.PD;
            subjects = dir(pd);

            for i = 1:length(subjects)
                subName = subjects(i).name;
                subf = fullfile(pd, subName);

                if ~isfolder(subf) || startsWith(subName, '.') || startsWith(subName,"..")
                    continue;
                end
                diary(sprintf('%s_output_log.txt',subName));

                flp = fullfile(subf, 'firstlevel');
                if ~isfolder(flp)
                    fprintf('%s doesn''t exist in %s\n', flp, subf);
                    continue;
                end

                % Load condition files
                matfiles = dir(fullfile(flp, '*condition_*.mat'));
                if isempty(matfiles)
                    fprintf('No condition files found in %s\n', flp);
                    continue;
                end

                matfile_name = matfiles(1).name;  % Assuming the first condition file
                matfile = load(fullfile(flp, matfile_name));

                % Get preprocessed fMRI files
                fl_files = dir(fullfile(flp, sprintf('*%s*.nii',uqIdentifier)));

                for run_index = 1:length(fl_files)
                    afs = fl_files(run_index).name;
                    preproc_fmri_file = fullfile(flp, afs);
                    pfiles = spm_select('ExtFPList',flp,afs,Inf);
                    cd(flp);

                    if ~exist(preproc_fmri_file, 'file')
                        fprintf('Preprocessed file %s not found\n', preproc_fmri_file);
                        continue;
                    end

                    % Specifying the main configuration
                    spm('defaults', 'FMRI');
                    matlabbatch = {};
                    matlabbatch{1}.spm.stats.fmri_spec.dir =cellstr(flp);
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = obj.Units;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = obj.RT;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

                    % Check scans input
                    disp('Scans input:');
                    disp(pfiles);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(pfiles);

                    % Create conditions in the batch
                    for condIdx = 1:NumConditions
                        cond_name = cond_name_list(condIdx);
                        has_pmod = ismember(cond_name, pmods);

                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).name = char(cond_name);
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).onset = matfile.onsets{condIdx};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).duration = matfile.durations{condIdx};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).tmod = 0;

                        if has_pmod
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).name = char(cond_name);
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).param = matfile.pmod(condIdx).param;
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).poly = matfile.pmod(condIdx).poly;
                            ll = sprintf('pmod script working: %s',char(cond_name));
                            disp(ll);
                        else
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                            lg = sprintf('non pmod script working: %s',char(cond_name));
                            disp(lg);
                        end
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).orth = 1;
                    end

                    % Additional SPM settings
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
                    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
                    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

                    % Correct dependencies
                    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
                        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
                    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

                    % Generate contrasts
                    pairs = nchoosek(cond_name_list, 2);  % Generate all pairs
                    consess_index = 1;  % Initialize index for consess

                    for pair_idx = 1:size(pairs, 1)
                        cond1 = pairs(pair_idx, 1);
                        cond2 = pairs(pair_idx, 2);

                        contrast_vector = arrayfun(@(cond) 1 * strcmp(cond, cond1) - 1 * strcmp(cond, cond2), cond_name_list);
                        reverse_contrast_vector = -contrast_vector;
                        disp('------------------');
                        disp(contrast_vector);
                        disp(reverse_contrast_vector);
                        disp('------------------');

                        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', ...
                            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.name = sprintf('%s - %s', char(cond1), char(cond2));
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.weights = contrast_vector;
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.sessrep = 'none';
                        consess_index = consess_index + 1;

                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.name = sprintf('%s - %s', char(cond2), char(cond1));
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.weights = reverse_contrast_vector;
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.sessrep = 'none';
                        consess_index = consess_index + 1;
                    end

                    matlabbatch{3}.spm.stats.con.delete = 0;

                    % Run the batch
                    try
                        spm_jobman('run', matlabbatch);
                    catch ME
                        disp('Error running spm_jobman:');
                        disp(ME.message);
                    end
                    save(['first_level', afs, '.mat'], 'matlabbatch');
                    diary off;
                end
            end
        end

        function firstRunDefined(obj, NumConditions, cond_name_list, pmods, contrast_names)
            % Method: firstrun
            %
            % Runs SPM batch processing for each subject's fMRI data, including model
            % specification, GLM estimation, and contrast generation. The function creates
            % and runs a batch using the provided number of conditions, condition names,
            % parametric modulation (pmod) details, and contrast names.
            %
            % **Usage:**
            % ```
            % obj.firstrun(NumConditions, cond_name_list, pmods, contrast_names)
            % ```
            %
            % **Parameters:**
            % - `NumConditions` (int): Number of conditions in the fMRI experiment.
            % - `cond_name_list` (cell array of char): List of condition names.
            % - `pmods` (cell array of char): List of conditions that have parametric modulation.
            % - `contrast_names` (cell array of char): List of contrast names to compute.
            %   Only forward and reverse contrasts for these names will be generated.
            %
            % **Returns:**
            % - None

            pd = obj.PD;
            subjects = dir(pd);

            for i = 1:length(subjects)
                subName = subjects(i).name;
                subf = fullfile(pd, subName);

                if ~isfolder(subf) || startsWith(subName, '.')
                    continue;
                end

                flp = fullfile(subf, 'firstlevel');
                if ~isfolder(flp)
                    fprintf('%s doesn''t exist in %s\n', flp, subf);
                    continue;
                end

                % Load condition files
                matfiles = dir(fullfile(flp, '*condition_*.mat'));
                if isempty(matfiles)
                    fprintf('No condition files found in %s\n', flp);
                    continue;
                end

                matfile_name = matfiles(1).name;  % Assuming the first condition file
                matfile = load(fullfile(flp, matfile_name));

                % Get preprocessed fMRI files
                fl_files = dir(fullfile(flp, '*swarADHD*.nii'));

                for run_index = 1:length(fl_files)
                    afs = fl_files(run_index).name;
                    preproc_fmri_file = fullfile(flp, afs);
                    pfiles = spm_select('ExtFPList', flp, afs, Inf);
                    cd(flp);

                    if ~exist(preproc_fmri_file, 'file')
                        fprintf('Preprocessed file %s not found\n', preproc_fmri_file);
                        continue;
                    end

                    % Specifying the main configuration
                    spm('defaults', 'FMRI');
                    matlabbatch = {};
                    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(flp);
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = obj.Units;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = obj.RT;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

                    % Check scans input
                    disp('Scans input:');
                    disp(pfiles);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(pfiles);

                    % Create conditions in the batch
                    for condIdx = 1:NumConditions
                        cond_name = cond_name_list(condIdx);
                        has_pmod = ismember(cond_name, pmods);

                        selected_durations = cell(1, length(condition_names));

                        % Loop through each condition name provided as input
                        for i = 1:length(condition_names)
                            % Find the index of the current condition name in the 'name' field of the struct
                            idx = find(strcmp(s.name, condition_names{i}));

                            % Check if the condition exists in the struct 'name' field
                            if ~isempty(idx)
                                % Get the corresponding duration using the index
                                selected_durations{i} = s.duration{idx};
                            else
                                % If condition is not found, display a warning and continue
                                warning('Condition "%s" not found in the struct "s".', condition_names{i});
                                selected_durations{i} = [];
                            end
                        end
                    end

                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).name = char(cond_name);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).onset = matfile.onsets{condIdx};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).duration = matfile.durations{condIdx};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).tmod = 0;

                    if has_pmod
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).name = char(cond_name);
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).param = matfile.pmod(condIdx).param;
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).poly = matfile.pmod(condIdx).poly;
                        ll = sprintf('pmod script working: %s', char(cond_name));
                        disp(ll);
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        lg = sprintf('non pmod script working: %s', char(cond_name));
                        disp(lg);
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).orth = 1;
                end

                % Additional SPM settings
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

                % Correct dependencies
                matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
                    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
                matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

                % Generate contrasts for specified contrast names
                consess_index = 1;  % Initialize index for consess
                for pair_idx = 1:length(contrast_names)
                    cond1 = contrast_names(pair_idx);
                    if ~ismember(cond1, cond_name_list)
                        warning('Contrast name "%s" not found in conditions. Skipping.', cond1);
                        continue;
                    end
                    remaining_conds = setdiff(cond_name_list, cond1);

                    for cond2_idx = 1:length(remaining_conds)
                        cond2 = remaining_conds(cond2_idx);

                        contrast_vector = arrayfun(@(cond) 1 * strcmp(cond, cond1) - 1 * strcmp(cond, cond2), cond_name_list);
                        reverse_contrast_vector = -contrast_vector;

                        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', ...
                            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.name = sprintf('%s - %s', cond1, cond2);
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.weights = contrast_vector;
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.sessrep = 'none';
                        consess_index = consess_index + 1;

                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.name = sprintf('%s - %s', cond2, cond1);
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.weights = reverse_contrast_vector;
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.sessrep = 'none';
                        consess_index = consess_index + 1;
                    end
                end

                matlabbatch{3}.spm.stats.con.delete = 0;

                % Run the batch
                try
                    spm_jobman('run', matlabbatch);
                catch ME
                    disp('Error running spm_jobman:');
                    disp(ME.message);
                end
                save(['first_level', afs, '.mat'], 'matlabbatch');
            end
        end
    end
end



%first level func from the reiter class - to be preserved
function firstRunAll(obj, NumConditions, cond_name_list, pmods,uqIdentifier)
            % Method: firstrun
            %
            % Runs SPM batch processing for each subject's fMRI data, including model
            % specification, GLM estimation, and contrast generation. The function creates
            % and runs a batch using the provided number of conditions, condition names,
            % and parametric modulation (pmod) details.
            %
            % **Usage:**
            % ```
            % obj.firstrun(NumConditions, cond_name_list, pmods)
            % ```
            %
            % **Parameters:**
            % - `NumConditions` (int): Number of conditions in the fMRI experiment.
            % - `cond_name_list` (cell array of char): List of condition names.
            % - `pmods` (cell array of char): List of conditions that have parametric modulation.
            %
            % **Returns:**
            % - None
            pd = obj.PD;
            subjects = dir(pd);

            for i = 1:length(subjects)
                subName = subjects(i).name;
                subf = fullfile(pd, subName);

                if ~isfolder(subf) || startsWith(subName, '.') || startsWith(subName,"..")
                    continue;
                end
                

                flp = fullfile(subf, 'firstlevel');
                if ~isfolder(flp)
                    fprintf('%s doesn''t exist in %s\n', flp, subf);
                    continue;
                end

                % Load condition files
                matfiles = dir(fullfile(flp, '*condition_*.mat'));
                if isempty(matfiles)
                    fprintf('No condition files found in %s\n', flp);
                    continue;
                end

                matfile_name = matfiles(1).name;  % Assuming the first condition file
                matfile = load(fullfile(flp, matfile_name));

                % Get preprocessed fMRI files
                fl_files = dir(fullfile(flp, sprintf('*%s*.nii',uqIdentifier)));

                for run_index = 1:length(fl_files)
                    afs = fl_files(run_index).name;
                    preproc_fmri_file = fullfile(flp, afs);
                    pfiles = spm_select('ExtFPList',flp,afs,Inf);
                    cd(flp);

                    diary(sprintf('%s_output_log.txt',subName));

                    if ~exist(preproc_fmri_file, 'file')
                        fprintf('Preprocessed file %s not found\n', preproc_fmri_file);
                        continue;
                    end

                    % Specifying the main configuration
                    spm('defaults', 'FMRI');
                    matlabbatch = {};
                    matlabbatch{1}.spm.stats.fmri_spec.dir =cellstr(flp);
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = obj.Units;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = obj.RT;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

                    % Check scans input
                    disp('Scans input:');
                    disp(pfiles);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(pfiles);

                    % Create conditions in the batch
                    for condIdx = 1:NumConditions
                        cond_name = cond_name_list(condIdx);
                        has_pmod = ismember(cond_name, pmods);

                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).name = char(cond_name);
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).onset = matfile.onsets{condIdx};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).duration = matfile.durations{condIdx};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).tmod = 0;

                        if has_pmod
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).name = char(cond_name);
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).param = matfile.pmod(condIdx).param;
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(1).poly = matfile.pmod(condIdx).poly;
                            ll = sprintf('pmod script working: %s',char(cond_name));
                            disp(ll);
                        else
                            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                            lg = sprintf('non pmod script working: %s',char(cond_name));
                            disp(lg);
                        end
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).orth = 1;
                    end

                    % Additional SPM settings
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
                    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
                    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

                    % Correct dependencies
                    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
                        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
                    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

                    % Generate contrasts
                    pairs = nchoosek(cond_name_list, 2);  % Generate all pairs
                    consess_index = 1;  % Initialize index for consess

                    for pair_idx = 1:size(pairs, 1)
                        cond1 = pairs(pair_idx, 1);
                        cond2 = pairs(pair_idx, 2);

                        contrast_vector = arrayfun(@(cond) 1 * strcmp(cond, cond1) - 1 * strcmp(cond, cond2), cond_name_list);
                        reverse_contrast_vector = -contrast_vector;
                        disp('------------------');
                        disp(contrast_vector);
                        disp(reverse_contrast_vector);
                        disp('------------------');

                        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', ...
                            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.name = sprintf('%s - %s', char(cond1), char(cond2));
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.weights = contrast_vector;
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.sessrep = 'none';
                        consess_index = consess_index + 1;

                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.name = sprintf('%s - %s', char(cond2), char(cond1));
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.weights = reverse_contrast_vector;
                        matlabbatch{3}.spm.stats.con.consess{consess_index}.tcon.sessrep = 'none';
                        consess_index = consess_index + 1;
                    end

                    matlabbatch{3}.spm.stats.con.delete = 0;

                    % Run the batch
                    try
                        spm_jobman('run', matlabbatch);
                    catch ME
                        disp('Error running spm_jobman:');
                        disp(ME.message);
                    end
                    save(['first_level', afs, '.mat'], 'matlabbatch');
                    diary off;
                end
            end
        end
