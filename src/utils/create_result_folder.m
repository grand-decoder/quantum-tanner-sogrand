function fullPath = create_result_folder(file_name_part)
% SAVETIMEDOBJECT Saves a variable to a timestamped folder and filename.
%
% Usage:
%   path = saveTimedObject(myVariable, 'myExperiment');
%   path = saveTimedObject(myVariable, 'myExperiment', 'C:\Data');
%
% Inputs:
%   inputObject - The data/object you want to save.
%   baseName    - The name for the file (and the variable inside the .mat).
    script_dir = fileparts(mfilename('fullpath'));
    root_folder = fullfile(script_dir, '../..');

    targetDir = fullfile(root_folder, "data", "output");

    % 2. Get current time once to ensure folder/file sync
    t = datetime('now');

    % 3. Create Folder Path (Format: yy_MM_dd)
    % Note: We use 'MM' for month. 'mm' is minutes in MATLAB.
    dateFolder = string(t, 'yy_MM_dd'); 
    saveDir = fullfile(targetDir, dateFolder);

    % Create directory if it doesn't exist
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    % 4. Create File Name (Format: HH-mm-ss)
    timeSuffix = string(t, 'HH_mm_ss');
    filename = timeSuffix + "_" + file_name_part;
    fullPath = fullfile(saveDir, filename);

    if ~exist(fullPath, 'dir')
        mkdir(fullPath);
    end
end