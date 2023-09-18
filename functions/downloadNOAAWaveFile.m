function outputFile = downloadNOAAWaveFile(url)
% Need to update output file name to match url name
    outputFile = 'gdaswave.t12z.gsouth.0p25.f000.grib2';  % Specify the name of the output file

    try
        % Download the file via HTTPS
        websave(outputFile, url);

        disp('File download successful.');
    catch
        % Handle errors
        disp('File download failed.');
    end
end