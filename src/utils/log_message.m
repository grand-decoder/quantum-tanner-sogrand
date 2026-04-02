function log_message(result_folder_path, message)
% Appends a timestamped message to log.txt in the specified folder.
%
% Arguments:
%   result_folder_path - Path to the folder where log.txt will be written
%   message            - Text string to record in the log

    log_file = fullfile(result_folder_path, "log.txt");
    try
        % Attempt to open file for appending
        file_id = fopen(log_file, 'a');
        % Check if file_id is valid (-1 means it failed to open)
        if file_id == -1
            warning('Could not open log file for writing: %s', log_file);
            return;
        end
        % Attempt to write
        time_stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        fprintf(file_id, '[%s] %s\n', time_stamp, message);
        % Close the handle
        fclose(file_id);
    catch me_exception
        % If something fails (e.g., disk full), MATLAB won't crash
        warning('Logging failed: %s', me_exception.message);
    end
end