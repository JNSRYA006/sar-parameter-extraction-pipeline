function fileSaveName = downloadNOAAWaveFile(url,fileSaveName)
    try
        % Download the file via HTTPS
        websave(fileSaveName, url);

        disp('File download successful.');
    catch
        % Handle errors
        disp('File download failed.');
    end
end