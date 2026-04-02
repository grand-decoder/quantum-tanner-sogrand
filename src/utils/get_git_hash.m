function git_hash = get_git_hash()
% Returns the current Git commit hash of the repository.
% Falls back to 'not available' if Git is not installed or the
% working directory is not a Git repository.
%
% Arguments:
%   git_hash - (output) Git commit hash string, or 'not available'

    % Default value if everything fails
    git_hash = 'not available';
    
    try
        % 1. Check if git even exists on the system path
        [status_check, ~] = system('git --version');
        if status_check ~= 0
            return; % Exit early with default value
        end
        
        % 2. Attempt to get the hash
        [status, cmdout] = system('git rev-parse HEAD');
        
        if status == 0
            git_hash = strtrim(cmdout);
        end
    catch
        % If system() itself throws an error, catch it here
        % git_hash remains 'not available'
    end
end