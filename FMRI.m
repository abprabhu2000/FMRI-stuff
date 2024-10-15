classdef FMRI
    properties
        PD           % Parent Directory
        NZS          % Number of Z Slices
        TR           % Repetition Time
        SO           % Slice Order
        gKernel      % Gaussian Smoothening Kernel
    end
    
    methods
        % Constructor
        function obj = FMRI(parentdir, numZSlices, repTime, sliceOrder,Gauss)
            obj.PD = parentdir;
            obj.NZS = numZSlices;
            obj.TR = repTime;
            obj.gKernel = Gauss;
            obj.SO = sliceOrder;
        end
        
        % Method to calculate repetition time (TA)
        function ta_vals = Reptimecalc(obj, funcfile)
            V = spm_vol(funcfile); % SPM function to get volume information
            zS = V(1).dim(3); % Number of slices in the z-dimension
            ta_vals = obj.TR - (obj.TR / zS); % Calculate TA
        end
        
        % Method to preprocess and execute SPM batch processing
        function Spm_PreprocExec(obj)
            results = containers.Map(); % Use a map to store results

            % Iterate over each subject directory
            subs = dir(obj.PD);
            for i = 1:length(subs)
                subName = subs(i).name;
                subject_path = fullfile(obj.PD, subName);

                % Skip hidden files/folders and non-directories
                if startsWith(subName, '.') || ~subs(i).isdir
                    continue;
                end

                % Initialize lists for T1 and func paths
                t1 = {};
                funcs = {};

                % Get the anat directory path
                anat_dir = fullfile(subject_path, 'anat');
                if isfolder(anat_dir)
                    t1_files = dir(fullfile(anat_dir, '*T1w.nii'));
                    if ~isempty(t1_files)
                        t1 = fullfile(anat_dir, t1_files(1).name);
                    end
                end
                disp(t1)
                % Get the func directory path
                func_dir = fullfile(subject_path, 'func');
                if isfolder(func_dir)
                    func_files = dir(fullfile(func_dir, '*bold.nii'));
                    for f = 1:length(func_files)
                        funcs{end+1} = fullfile(func_dir, func_files(f).name);
                    end
                end
                disp(funcs)
               
                if ischar(t1)
                    t1 = {t1};
                elseif ~iscell(t1)
                    t1 = cellstr(t1);
                end

                if ischar(funcs)
                    funcs = {funcs};
                elseif ~iscell(funcs)
                    funcs = cellstr(funcs);
                end

                % Ensure both t1 and funcs are column cell arrays
                t1 = t1(:);
                funcs = funcs(:);

                % Combine anat and func paths into a single cell array
                all_paths = [t1; funcs];
                
                disp('------------------');
                disp(all_paths);
                disp('------------------');


                % Add results to the map
                if ~isempty(all_paths)
                    results(subName) = all_paths;
                else
                    disp(['No files found for subject ', subName]);
                end
            end

            % Iterate over each subject and process the data
            keys = results.keys;
            for i = 1:length(keys)
                subs = keys{i};
                paths = results(subs);
                anat_p = paths{1}; % Assuming the first path is the anatomical image
                func_p = paths(2:end); % The rest are functional images

                for j = 1:length(func_p)
                    if exist(func_p{j}, 'file')
                        cd(fileparts(func_p{j}));

                        func_files = spm_select('ExtFPList', fileparts(func_p{j}), ['^', spm_file(func_p{j}, 'basename'), '.nii$'], Inf);
                        anat_files = spm_select('FPList', fileparts(anat_p), ['^', spm_file(anat_p, 'basename'), '.nii$']);
                        
                        % Prepare batch
                        matlabbatch = {};
                        matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cellstr(func_files);
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.95;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2; %<-- can be user input  / wait for lorenz's input (4/5)
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
                        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

                        matlabbatch{2}.spm.temporal.st.scans{1}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', ...
                            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
                        matlabbatch{2}.spm.temporal.st.nslices = obj.NZS;
                        matlabbatch{2}.spm.temporal.st.tr = obj.TR;
                        matlabbatch{2}.spm.temporal.st.ta = obj.Reptimecalc(func_p{j});
                        matlabbatch{2}.spm.temporal.st.so = obj.SO;
                        matlabbatch{2}.spm.temporal.st.refslice = 1; % can be aligned to mean / can refer to middle slice 
                        matlabbatch{2}.spm.temporal.st.prefix = 'a';

                        matlabbatch{3}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', ...
                            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
                        matlabbatch{3}.spm.spatial.coreg.estwrite.source = cellstr(anat_files);
                        matlabbatch{3}.spm.spatial.coreg.estwrite.other = {''};
                        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                        matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                        matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

                        % Check if the preprocessed files already exist
                        tpm_files_exist = exist('c1sub-01_T1w.nii', 'file') && ...
                            exist('c2sub-01_T1w.nii', 'file') && ...
                            exist('c3sub-01_T1w.nii', 'file') && ...
                            exist('c4sub-01_T1w.nii', 'file') && ...
                            exist('c5sub-01_T1w.nii', 'file') && ...
                            exist('msub-01_T1w.nii', 'file') && ...
                            exist('rsub-01_T1w.nii', 'file') && ...
                            exist('y_sub-01_T1w.nii', 'file');

                        if ~tpm_files_exist
                            matlabbatch{4}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate & Reslice: Coregistered Images', ...
                                substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
                            matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
                            matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
                            matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];
                            matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,1')};
                            matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
                            matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
                            matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
                            matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
                            matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
                            matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
                            matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
                            matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
                            matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
                            matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
                            matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
                            matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
                            matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
                            matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
                            matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
                            matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
                            matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
                            matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
                            matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
                            matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];
                            matlabbatch{4}.spm.spatial.preproc.warp.vox = NaN;
                            matlabbatch{4}.spm.spatial.preproc.warp.bb = [NaN NaN NaN; NaN NaN NaN];
                        else
                            disp('Segmented Structural files already exist, skipping spatial preprocessing step.'); 
                        end

                        % Continue with the remaining steps in the batch
                        matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', ...
                            substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
                        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', ...
                            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
                        matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
                        matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2]; % can be changed - resampled voxel size
                        matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
                        matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
                        matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', ...
                            substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
                        matlabbatch{6}.spm.spatial.smooth.fwhm = obj.gKernel; % fwhm smoothing kernel 2 - 4 times vox size 
                        matlabbatch{6}.spm.spatial.smooth.dtype = 0;
                        matlabbatch{6}.spm.spatial.smooth.im = 0;
                        matlabbatch{6}.spm.spatial.smooth.prefix = 's'; % use s6 , s8 as gaussian kernel prefix 

                        % Run the batch
                        spm_jobman('run', matlabbatch);

                        % Save the batch to a .mat file
                        save(['preprocessing_batch_' subs '_run_' num2str(j) '.mat'], 'matlabbatch');
                    end
                end
            end
        end
    end
end
