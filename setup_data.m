% Downloads quantum Tanner codes from the qTanner repository.

% Configuration
commitHash = '5966d702702610b8be98c50fc0443d1a54e98f96';
url = sprintf('https://github.com/RebKatRad/qTanner/archive/%s.zip', commitHash);
tempDir = 'temp_download';
finalDest = fullfile('data', 'codes');

if ~exist(finalDest, 'dir')
    mkdir(finalDest);
end

unzip(url, tempDir);
srcPath = fullfile(tempDir, ['qTanner-', commitHash], 'data');

% move code folder to data directory
movefile(fullfile(srcPath, '*'), finalDest);

% Cleanup the temporary directory
rmdir(tempDir, 's');
