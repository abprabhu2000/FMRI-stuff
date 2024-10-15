classdef firstLev2
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
    % 3. **pmod generator:** Creates PMOD files from the SPM condition files.
    % 4. **firstrun:** Runs the SPM batch processing for each subject's fMRI data.
    % 5. **FirstLevelCleaner:** Cleans up the output files by organizing them into
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
        function obj = firstLevReiter(parent, RT, units)
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
            % **Parameters:*
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
                        if ~isempty(regexp(fileName,'^rp_sub-.*_bold.txt$','once')) && endsWith(filename,'.txt')
                            continue;
                        else
                            req_file_folder = fullfile(subf, 'func');
                            if isfolder(req_file_folder)
                                funcs = dir(req_file_folder);
                                for k = 1:length(funcs)
                                    funcName = funcs(k).name;
                                    if ~isempty(regexp(funcName,'^rp_sub-.*_bold.txt$','once')) && endsWith(funcName,'.txt')
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

        function [pmods_cat] = pmodCatGen(obj) %#ok<MANU>
            % PMODCATGEN generates a struct with categorized names for given parametric modulations.
            %
            % This method searches for '*condition_*.mat' files in the specified folder path (flp),
            % loads the parametric modulation values, categorizes them, and saves the results in a new .mat file.
            %
            % Parameters:
            %   - obj: The class object.
            %   - flp: Folder path where condition .mat files are located.
            %
            % Returns:
            %   - pmods_cat: A struct with categorized fields ('win', 'loss', 'late') and corresponding values.

            % Navigate to the specified folder
            % Find all condition .mat files in the directory
            matfiles = dir('*condition_*.mat');

            % Check if valid files are found (exclude files that start with '.' or '..')
            matfile_selected = [];
            for i = 1:length(matfiles)
                if ~startsWith(matfiles(i).name, '.')
                    matfile_selected = matfiles(i);  % Select the first valid file
                    break;  % Stop after finding the first valid file
                end
            end

            if isempty(matfile_selected)
                fprintf('No valid condition files found in %s\n', flp);
                pmods_cat = [];  % Return empty if no valid file is found
                return;
            end

            % Load the selected condition file
            matfile = load(matfile_selected.name);

            % Assuming the 4th index has the relevant pmod
            if isfield(matfile, 'pmod') && length(matfile.pmod) >= 4 && isfield(matfile.pmod(4), 'param')
                pmod_arr = matfile.pmod(4).param;
            else
                fprintf('The selected file does not have the expected pmod structure.\n');
                pmods_cat = [];
                return;
            end

            % Initialize categorized arrays
            pwin = [];
            ploss = [];
            plate = [];

            % Loop through each row in the input matrix and categorize based on the values in the array
            for i = 1:size(pmod_arr, 1)
                p_arr = pmod_arr(i, :);  % Extract the current row

                % Categorize based on the values in the array
                pwin(i, :) = p_arr .* (p_arr == 1);    % Extract values equal to 1
                ploss(i, :) = p_arr .* (p_arr == -1);  % Extract values equal to -1
                plate(i, :) = p_arr .* (p_arr == 0.1); % Extract values equal to 0.1
            end

            % Create a struct to hold the categorized data
            pmods_cat = struct('conds', {'win', 'loss', 'late'}, 'vals', {pwin, ploss, plate});

            % Save the categorized data in a .mat file in the 'firstlevel' folder
            save(fullfile(pwd, "pmod_sep.mat"), "pmods_cat");
        end


        function pmodCatCreator(obj)
            % Method to create parametric modulation categories for all subjects.
            %
            % **Usage:**
            % ```
            % obj.pmodCatCreator()
            % ```
            %
            % **Parameters:**
            % - None

            pd = obj.PD;
            subjects = dir(pd);

            for i = 1:length(subjects)
                subName = subjects(i).name;
                subf = fullfile(pd, subName);

                if ~isfolder(subf) || startsWith(subName, '.') || startsWith(subName, "..")
                    continue;
                end

                flp = fullfile(subf, 'firstlevel');
                if ~isfolder(flp)
                    fprintf('%s doesn''t exist in %s\n', flp, subf);
                    continue;
                else
                    cd(flp);  % Change to the 'firstlevel' folder
                    obj.pmodCatGen();  % Call the function with the folder path
                end
            end
        end



        function firstRunModded(obj, NumConditions, cond_name_list, pmods,uqIdentifier)
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
            % - `uqIdentifier` : string that contains somethig that is unique to the preprocessed bold files like 'swarADHD' for example.
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
                
                % load motion params
                allFiles = dir(flp);
                motion_params = allFiles(~cellfun('isempty',regexp({allFiles.name},'^rp_sub-.*_bold.txt$')));
                matchedFiles = fullfile(flp,motion_params(1).name);

                % Load condition files
                matfiles_pmods = dir(fullfile(flp, 'pmod_sep.mat'));
                if isempty(matfiles_pmods)
                    fprintf('No pmod files found in %s\n', flp);
                    continue;
                end

                pmods_file = load(fullfile(flp,matfiles_pmods(1).name));
                pmvals = pmods_file.pmods_cat;

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
                    diaryFile = sprintf('%s_output_log.txt', subName);  % Generate the diary file name
                    if exist(diaryFile, 'file')  % Check if the file exists
                        delete(diaryFile);  % Delete the existing diary file
                    end
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
                            for idx = 1:length(pmvals)
                                pmvals(idx).vals(pmvals(idx).vals == 0) = abs(pmvals(idx).vals(pmvals(idx).vals == 0));
                                fprintf('\n-------------------- %s -------------------\n',pmvals(idx).conds);
                                fprintf('[%f] ',unique(pmvals(idx).vals));
                                fprintf('\n--------------------end--------------------\n');


                                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(idx).name = sprintf('%s',pmvals(idx).conds);
                                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(idx).param = pmvals(idx).vals;
                                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod(idx).poly = 1;
                                ll = sprintf('pmod script working: %s , pmod_name = %s',char(cond_name),pmvals(idx).conds);
                                disp(ll);
                            end
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
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(matchedFiles);
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
        function feedbackPrep(obj)
            % feedbackPrep - Prepares feedback interaction regressors for fMRI analysis.
            %
            % This function searches for each subject's SPM.mat file in their 'firstlevel' folder,
            % extracts the feedback interaction regressors from the design matrix, and saves them
            % as a structure ('fbStruct.mat') containing the names and values of the regressors.
            %
            % **Usage:**
            % ```
            % obj.feedbackPrep()
            % ```
            %
            % **Parameters:**
            % - `obj` (FirstLevel object): The instance of the class containing the subject data.
            %
            % **Behavior:**
            % - For each subject:
            %   1. Checks if the subject's 'firstlevel' folder exists.
            %   2. Searches for the 'SPM.mat' file in the 'firstlevel' folder.
            %   3. If found, loads the SPM.mat file and extracts feedback interaction regressors
            %      from the design matrix (xX.X).
            %   4. Creates a structure `fbStruct` with fields for the feedback interaction regressors.
            %   5. Saves the `fbStruct` in the 'firstlevel' folder as 'fbStruct.mat'.
            %
            % **Returns:**
            % - None
            %
            % **Example:**
            % ```
            % obj.feedbackPrep();
            % ```
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

                cd(flp)
                spm_mat = dir(fullfile(flp,'SPM.mat'));
                
                if ~isempty(spm_mat)
                    disp('SPM mat found');
                else
                    disp('SPM mat not found');
                    continue;
                end 
                spmM = load(spm_mat.name);
                spmM = spmM.SPM;
                fbStruct = struct('conds',{'feedbackxwin','feedbackxloss','feedbackxlate'},'vals',{spmM.xX.X(:,5),spmM.xX.X(:,6),spmM.xX.X(:,7)});
                save("fbStruct.mat",'fbStruct');
            end

        end

        function fbIntRun(obj,uqIdentifier)
            % fbIntRun - Runs SPM analysis on feedback interaction regressors for each subject.
            %
            % This function performs fMRI analysis by running SPM batch processing on feedback
            % interaction regressors for each subject. It loads feedback interaction data from the
            % saved 'fbStruct.mat' file and runs model specification, GLM estimation, and contrast
            % generation.
            %
            % **Usage:**
            % ```
            % obj.fbIntRun()
            % ```
            %
            % **Parameters:**
            % - `obj` (FirstLevel object): The instance of the class containing the subject data.
            %
            % **Behavior:**
            % - For each subject:
            %   1. Checks if the 'feedback' directory exists within the 'firstlevel' folder.
            %   2. Creates the 'feedback' directory if it does not exist.
            %   3. Loads the 'fbStruct.mat' file, which contains the feedback interaction data.
            %   4. Searches for the realignment parameters file (e.g., 'rp_sub-*_bold.txt') in the 'firstlevel' folder.
            %   5. Loads preprocessed fMRI files based on the unique identifier.
            %   6. Specifies the SPM batch configuration, including model specification, GLM estimation,
            %      and contrast generation.
            %   7. Generates all possible pairwise contrasts for feedback interaction regressors.
            %   8. Runs the SPM batch using `spm_jobman`.
            %   9. Saves the batch configuration for each subject.
            %
            % **Returns:**
            % - None
            %
            % **Example:**
            % ```
            % obj.fbIntRun();
            % ```
            pd = obj.PD;
            subjects = dir(pd);

            for i = 1:length(subjects)
                subName = subjects(i).name;
                subf = fullfile(pd, subName);

                if ~isfolder(subf) || startsWith(subName, '.') || startsWith(subName,"..")
                    continue;
                end


                flp = fullfile(subf, 'firstlevel');
                

                if ~isfolder('feedback')
                    mkdir(fullfile(flp,'feedback'));
                else
                    disp('feedback dir already exists');
                end 

                if ~isfolder(fullfile(flp,'feedback'))
                    fprintf('%s doesn''t exist in %s\n', flp, subf);
                    continue;
                end

                % load fbstruct
                
                if isfolder(fullfile(flp,'feedback'))
                    fb_filepath = fullfile(flp,'feedback');
                else
                    disp('feedback folder not found');
                end
                feedback = dir(fullfile(flp,'fbstruct.mat'));
                fb = load(fullfile(feedback.folder,feedback.name));
                fb = fb.fbStruct;
                
                allFiles = dir(flp);
                motion_params = allFiles(~cellfun('isempty',regexp({allFiles.name},'^rp_sub-.*_bold.txt$')));
                matchedFiles = fullfile(flp,motion_params(1).name);
                

                % Get preprocessed fMRI files
                fl_files = dir(fullfile(flp, sprintf('*%s*.nii',uqIdentifier)));

                for run_index = 1:length(fl_files)
                    afs = fl_files(run_index).name;
                    preproc_fmri_file = fullfile(flp, afs);
                    pfiles = spm_select('ExtFPList',flp,afs,Inf);
                    cd(fb_filepath);
                    diaryFile = sprintf('%s_output_log_feedback.txt', subName);  % Generate the diary file name
                    if exist(diaryFile, 'file')  % Check if the file exists
                        delete(diaryFile);  % Delete the existing diary file
                    end
                    diary(sprintf('%s_output_log_feedback.txt',subName));

                    if ~exist(preproc_fmri_file, 'file')
                        fprintf('Preprocessed file %s not found\n', preproc_fmri_file);
                        continue;
                    end

                    % Specifying the main configuration
                    spm('defaults', 'FMRI');
                    matlabbatch = {};
                    matlabbatch{1}.spm.stats.fmri_spec.dir =cellstr(fb_filepath);
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = obj.Units;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = obj.RT;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

                    % Check scans input
                    disp('Scans input:');
                    disp(pfiles);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(pfiles);

                    % Create conditions in the batch
                    for condIdx = 1:length(fb)
                        cond_name = fb(condIdx).conds;
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).name = char(cond_name);
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).onset = fb(condIdx).vals;
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).duration = 0;
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).tmod = 0;


                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        lg = sprintf('script working: %s',char(cond_name));
                        disp(lg);

                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(condIdx).orth = 1;
                    end

                    % Additional SPM settings
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(matchedFiles);
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
                    pairs = nchoosek({fb.conds}, 2);  % Generate all pairs
                    consess_index = 1;  % Initialize index for consess

                    for pair_idx = 1:size(pairs, 1)
                        cond1 = pairs(pair_idx, 1);
                        cond2 = pairs(pair_idx, 2);

                        contrast_vector = arrayfun(@(cond) 1 * strcmp(cond, cond1) - 1 * strcmp(cond, cond2), {fb.conds});
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

    end
end

