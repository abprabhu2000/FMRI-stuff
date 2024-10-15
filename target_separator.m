function process_all_files(directory)
    % Get a list of all .mat files in the directory
    matFiles = dir(fullfile(directory, '*.mat'));
    
    % Loop through each .mat file in the directory
    for i = 1:length(matFiles)
        matFile = fullfile(directory, matFiles(i).name);
        fprintf('Processing %s...\n', matFile);
        process_mat_file(matFile);
    end
end

function process_mat_file(mat_file)
    % Load the .mat file
    matData = load(mat_file);
    onsets = matData.onsets;
    
    % Extract cue, target, error, and feedback onsets
    cue = onsets{1}(:)';
    target = onsets{2}(:)';
    error = onsets{3}(:)';
    feedback = onsets{4}(:)';
    
    % Step 1: Perform array math, concatenate error and feedback
    feed_full = sort([error feedback]);
    
    % Step 2: Insert NaNs from index 20 to 40
    feed_full_new = insert_nans(feed_full, 20, 40);

    % Step 3: Calculate differences
    [cfb, tfb] = IIT_calc(cue, target, feed_full_new);
    
    % Step 4: Save the results to text files
    save_differences(mat_file, cfb, tfb);
end

function new_arr = insert_nans(arr, start_idx, end_idx)
    % Create a new array with space for NaNs
    new_length = length(arr) + (end_idx - start_idx);
    new_arr = nan(1, new_length);
    
    % Insert values from the original array before the 'start_idx'
    new_arr(1:start_idx) = arr(1:start_idx);
    
    % Insert values after 'start_idx' into the new array starting from 'end_idx'
    new_arr(end_idx+1:end) = arr(start_idx+1:end);
end

function [result_cfb, result_tfb] = IIT_calc(cue, target, feedback)
    % Remove the first cue and target, and the last feedback
    nc = cue(2:end);
    nt = target(2:end);
    nfb = feedback(1:end-1);
    
    result_cfb = nan(1, length(nc));
    result_tfb = nan(1, length(nt));
    
    % Iterate over the valid feedback and cue values
    for i = 1:length(result_cfb)
        if ~isnan(nfb(i))  % Ensure the feedback is not NaN
            if ~isnan(nc(i))
                result_cfb(i) = nc(i) - nfb(i);  % Difference between next cue and last feedback
            end
            if ~isnan(nt(i))
                result_tfb(i) = nt(i) - nfb(i);  % Difference between next target and last feedback
            end
        end
    end
end

function save_differences(mat_file, cfb, tfb)
    [~, base_name, ~] = fileparts(mat_file);  % Get the base name of the file
    
    % Save Feedback-Cue differences
    cfb_file = fullfile(pwd, [base_name '_cfb.txt']);
    fid = fopen(cfb_file, 'w');
    fprintf(fid, 'Feedback Cue Diff\n');
    fprintf(fid, '%f\n', cfb);
    fclose(fid);
    
    % Save Feedback-Target differences
    tfb_file = fullfile(pwd, [base_name '_tfb.txt']);
    fid = fopen(tfb_file, 'w');
    fprintf(fid, 'Feedback Target Diff\n');
    fprintf(fid, '%f\n', tfb);
    fclose(fid);
end

function [target_error, target_feedback] = process_feedback_error_target(feedback, error, target)
    % Create error and feedback masks
    error_mask = [error, ones(length(error), 1)];
    feedback_mask = [feedback, zeros(length(feedback), 1)];
    
    % Concatenate and sort by the first column
    tot_feedback = [error_mask; feedback_mask];
    [~, idx] = sort(tot_feedback(:, 1));
    sort_tot = tot_feedback(idx, :);
    
    % Extract the mask array
    mask_arr = sort_tot(:, 2);
    
    % Insert NaNs between index 20 and 40
    mask_arr = insert_nans(mask_arr, 20, 40);
    
    % Calculate target_error
    target_error = calculate_target_array(target, mask_arr);
    
    % Invert mask_arr (0 becomes 1, 1 becomes 0)
    mask_arr_inv = mask_arr;
    mask_arr_inv(mask_arr == 0) = 1;
    mask_arr_inv(mask_arr == 1) = 0;
    
    % Calculate target_feedback
    target_feedback = calculate_target_array(target, mask_arr_inv);
    
    % Remove NaNs from the results
    target_error = target_error(~isnan(target_error));
    target_feedback = target_feedback(~isnan(target_feedback));
end

function result = calculate_target_array(target, mask_arr)
    % Calculate the target array (target_error or target_feedback) based on the mask array
    result = target .* mask_arr;
    result(mask_arr == 0 | isnan(mask_arr)) = NaN; % Set invalid values to NaN
    result = result(~isnan(result)); % Remove NaNs
end

