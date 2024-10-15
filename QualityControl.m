classdef QualityControl
    properties
        PD      % Parent directory
        uqIA    % unique Identifier Anat
        uqIF    % unique identifier Func
    end

    methods
        % Constructor
        function obj = QualityControl(ParentDirectory, uniqueIdAnat, uniqueIdFunc)
            obj.PD = ParentDirectory;
            obj.uqIA = uniqueIdAnat;
            obj.uqIF = uniqueIdFunc;
        end

        function infoStruct = qcInfo(obj)
            % Initialize infoStruct with anatPaths and funcPaths as empty cell arrays
            infoStruct = struct('ParentDirectory', {obj.PD}, ...
                'uniqueIdAnat', {obj.uqIA}, ...
                'UniqueIdFunc', {obj.uqIF}, ...
                'anatPaths', {{}}, ...
                'funcPaths', {{}});
        end

        function [infoStruct, obj] = qcPreprocessing(obj, infoStruct)
            pd = obj.PD;
            subs = dir(pd);

            for i = 1:length(subs)
                sName = subs(i).name;
                if startsWith(sName, '.') || startsWith(sName, '_') || ~isfolder(fullfile(subs(i).folder, sName))
                    fprintf('%s not in subs folder, will not be included in the QC check \n', sName);
                    continue;
                end
                sublets = dir(fullfile(subs(i).folder, sName));

                % Look for anat subdirectory and process
                for anats = 1:length(sublets)
                    sublName = sublets(anats).name;
                    subletDir = fullfile(subs(i).folder, sName, sublName);
                    if startsWith(sublName, '.') || startsWith(sublName, '_') || ~isfolder(subletDir)
                        fprintf('%s in %s folder will be skipped \n', sublName, sName);
                        continue;
                    end
                    if strcmpi(sublName, 'anat')
                        allAnat = dir(fullfile(subletDir, '*.nii'));
                        fileNames = {allAnat.name};
                        uID = sprintf('^%s.*\\_T1w.nii$', obj.uqIA);
                        matchingFiles = fileNames(~cellfun('isempty', regexp(fileNames, uID)));

                        % Add matching files to anatPaths in infoStruct
                        for fileIdx = 1:length(matchingFiles)
                            fullFilePath = fullfile(subletDir, matchingFiles{fileIdx});  % Full path to the file
                            infoStruct.anatPaths{end+1} = fullFilePath;  % Append full path to anatPaths
                            infoStruct.anatPaths = infoStruct.anatPaths';
                        end

                    end
                end

                % Look for func subdirectory and process
                for funcs = 1:length(sublets)
                    sublName = sublets(funcs).name;
                    subletDir = fullfile(subs(i).folder, sName, sublName);
                    if startsWith(sublName, '.') || startsWith(sublName, '_') || ~isfolder(subletDir)
                        fprintf('%s in %s folder will be skipped \n', sublName, sName);
                        continue;
                    end
                    if strcmpi(sublName, 'func')
                        allFunc = dir(fullfile(subletDir, '*.nii'));
                        fileNames = {allFunc.name};
                        uID = sprintf('^%s.*\\_bold.nii$', obj.uqIF);
                        matchingFiles = fileNames(~cellfun('isempty', regexp(fileNames, uID)));

                        % Add matching files to funcPaths in infoStruct
                        for fileIdx = 1:length(matchingFiles)
                            fullFilePath = fullfile(subletDir, matchingFiles{fileIdx});  % Full path to the file
                            infoStruct.funcPaths{end+1} = fullFilePath;
                            infoStruct.funcPaths = infoStruct.funcPaths';
                        end

                    end
                end
            end
        end

        function qcRun(obj, infoStruct)
            infs = infoStruct;

            if length(infs.anatPaths) == length(infs.funcPaths)
                for idx = 1:length(infs.anatPaths)
                    % Extract subject name from anatomical image path
                    anatPath = infs.anatPaths{idx};  % Anatomical image path
                    funcPath = infs.funcPaths{idx};  % Functional image path

                    % Use regexp to extract the part matching 'sub-XX' where XX is any number
                    pattern = 'sub-\d{2}';  % Pattern to match 'sub-' followed by exactly two digits
                    match = regexp(anatPath, pattern, 'match');  % Find the match

                    if ~isempty(match)
                        subjectName = match{1};  % Extracted subject part, e.g., 'sub-03'
                    else
                        subjectName = 'Unknown Subject';  % Fallback in case of no match
                    end

                    % Prepare the images to display
                    imgs = char(anatPath, [funcPath ',1']);  % Only first volume of functional image

                    % Display a loading message
                    fprintf('Loading NIfTI file for %s , %s\n', anatPath, funcPath);

                    % Display both images in SPM
                    spm_check_registration(imgs);

                    % Add captions using the extracted subject name
                    spm_orthviews('Caption', 1, sprintf('%s: Anat Image', subjectName));
                    spm_orthviews('Caption', 2, sprintf('%s: First Func Image', subjectName));

                    % Wait for the user to close the SPM image window
                    waitfor(spm_figure('FindWin', 'Graphics'));  %<--------- unchanged
                end
            else
                fprintf('Mismatch in anat and func files.\n');
            end
        end


        function [infoStruct, obj] = qcRun2(obj, infoStruct)
            infs = infoStruct;  % Retrieve the infoStruct passed to the function

            % Initialize a new struct field for storing registration feedback
            regFeedback = struct('Subject', {}, 'QualCheck', {});

            if length(infs.anatPaths) == length(infs.funcPaths)
                for idx = 1:length(infs.anatPaths)
                    % Extract subject name from anatomical image path
                    anatPath = infs.anatPaths{idx};  % Anatomical image path
                    funcPath = infs.funcPaths{idx};  % Functional image path

                    % Use regexp to extract the part matching 'sub-XX' where XX is any number
                    pattern = 'sub-\d{2}';  % Pattern to match 'sub-' followed by exactly two digits
                    match = regexp(anatPath, pattern, 'match');  % Find the match

                    if ~isempty(match)
                        subjectName = match{1};  % Extracted subject part, e.g., 'sub-03'
                    else
                        subjectName = 'Unknown Subject';  % Fallback in case of no match
                    end

                    % Prepare the images to display
                    imgs = char(anatPath, [funcPath ',1']);  % Only first volume of functional image

                    % Display a loading message
                    fprintf('Loading NIfTI file for %s , %s\n', anatPath, funcPath);

                    % Display both images in SPM
                    spm_check_registration(imgs);

                    % Add captions using the extracted subject name
                    spm_orthviews('Caption', 1, sprintf('%s: Anat Image', subjectName));
                    spm_orthviews('Caption', 2, sprintf('%s: First Func Image', subjectName));

                    % Wait for the user to close the SPM image window
                    waitfor(spm_figure('FindWin', 'Graphics'));  %<--------- unchanged

                    % Prompt the user to confirm if they are satisfied with the registration
                    userResponse = questdlg(sprintf('Are you satisfied with the registration for %s?', subjectName), ...
                        'Registration Check', 'Yes', 'No', 'Yes');

                    % Record the user response
                    regFeedback(end+1).Subject = subjectName;
                    regFeedback(end).QualCheck = userResponse;
                end

                % Add the registration feedback to the infoStruct
                infoStruct.registrationFeedback = regFeedback;

                % Convert the registrationFeedback struct to a table
                feedbackTable = struct2table(infoStruct.registrationFeedback);

                % Define the filename for the Excel file
                excelFileName = fullfile(obj.PD,'pre-preprocessing QC check.xlsx');

                % Write the table to an Excel file
                writetable(feedbackTable, excelFileName);

                fprintf('Registration feedback saved to %s\n', excelFileName);
            else
                fprintf('Mismatch in anat and func files.\n');
            end
        end

        function process_movement_params(obj)
            parent_dir = obj.PD;
            % This function reads rp_*.txt files, computes basic stats, plots movement
            % parameters, and saves results to an Excel file and image.

            % Go through all subject directories
            subjects = dir(parent_dir);
            for i = 1:length(subjects)
                if subjects(i).isdir && ~startsWith(subjects(i).name, '.') && ~startsWith(subjects(i).name, '_')
                    subject_dir = fullfile(parent_dir, subjects(i).name);

                    % Find the 'func' directory within the subject's folder
                    func_dir = fullfile(subject_dir, 'func');
                    if isfolder(func_dir)
                        txt_files = dir(fullfile(func_dir, 'rp_*.txt')); % Find the rp_*.txt file

                        if ~isempty(txt_files)
                            txt_file = fullfile(func_dir, txt_files(1).name);

                            % Read the data from the text file
                            data = readmatrix(txt_file);

                            % Set the column names
                            columns = {'X', 'Y', 'Z', 'Roll', 'Yaw', 'Pitch'};

                            % Compute basic statistics
                            gg = struct();
                            for j = 1:6
                                gg.mean.(columns{j}) = mean(data(:, j));
                                gg.min.(columns{j}) = min(data(:, j));
                                gg.max.(columns{j}) = max(data(:, j));
                            end

                            % Compute ranges
                            ranges = struct();
                            for j = 1:6
                                ranges.(columns{j}) = gg.max.(columns{j}) - gg.min.(columns{j});
                            end

                            % Create a new table for Excel with the stats and ranges
                            T = struct2table(gg);
                            range_table = struct2table(ranges);
                            T = [T; range_table];

                            % Write the table to Excel
                            [~, name, ~] = fileparts(txt_file);
                            excel_filename = fullfile(func_dir, [name '.xlsx']);
                            writetable(T, excel_filename);
                            disp(['Excel file written: ', excel_filename]);

                            % Create the plots
                            figure('Position', [100, 100, 1000, 1000]);
                            subplot_titles = {'X movement param', 'Y movement param', 'Z movement param', ...
                                'Roll movement param', 'Yaw movement param', 'Pitch movement param'};

                            for j = 1:6
                                subplot(3, 2, j);
                                plot(1:size(data, 1), data(:, j));
                                hold on;
                                title(subplot_titles{j});
                                yline(gg.mean.(columns{j}), '--k', 'mean');
                                yline(gg.max.(columns{j}), '--r', 'max');
                                yline(gg.min.(columns{j}), '--g', 'min');
                                legend('Data', 'Mean', 'Max', 'Min');
                                hold off;
                            end

                            % Save the plot
                            image_filename = fullfile(func_dir, [name '.jpeg']);
                            saveas(gcf, image_filename);
                            disp(['Image file written: ', image_filename]);
                            close(gcf);
                        end
                    end
                end
            end
        end



    end
end
