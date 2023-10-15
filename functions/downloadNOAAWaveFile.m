function fileSaveName = downloadNOAAWaveFile(url,fileSaveName)
% Need to update output file name to match url name
    % For wave data
    %outputFile = 'gdaswave.t12z.gsouth.0p25.f000.grib2';  % Specify the name of the output file 
    % For sla
    %outputFile = 'download.nc';
    %options =
    %weboptions('Username','ryanjonesza','Password','vAr$N2PWW#@q#rm'); For
    %NASA
    try
        % Download the file via HTTPS
        websave(fileSaveName, url);

        disp('File download successful.');
    catch
        % Handle errors
        disp('File download failed.');
    end
end