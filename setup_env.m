% setup python environment.
try
    % Get the root folder where this script lives
    [projectRoot, ~, ~] = fileparts(mfilename('fullpath'));

    envPath = fullfile(projectRoot, '.venv', 'Scripts', 'python.exe');
    
    fprintf('Setting up Python environment: %s\n', envPath);
    
    % Check if the venv actually exists before trying to load it
    if exist(envPath, 'file')
        pyenv('Version', envPath);
    else
        error('Python executable not found at: %s. Did you create the venv on this machine?', envPath);
    end

catch ME
    warning('Failed to set up Python environment automatically.');
    fprintf('Error: %s\n', ME.message);
end

%
addpath("experiments");
addpath("generate_paper_plots");
addpath("src");
addpath("src/decoder");
addpath("src/utils");