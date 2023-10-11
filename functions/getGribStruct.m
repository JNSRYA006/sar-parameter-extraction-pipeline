function struc = getGribStruct(wgrib2Path,GribPath)

OutputColumns = ["latitude", "longitude", "time", "significantWaveHeight", "significantWavePeriod", "direction", "windSpeed", "windDirection"];
str = dir(strcat(GribPath,'\*.grib2'));
S = struct();
for i = 1:length(OutputColumns)
    S = setfield(S,OutputColumns(i),[]);
end

for i = 1:length(str)
    % Define the directory where wgrib2 is located
    wgrib2_dir = wgrib2Path;
    
    % Define the full path to the wgrib2 executable
    wgrib2_executable = fullfile(wgrib2_dir, 'wgrib2.exe');
    
    % Define the command to run
    file_path = strcat(str(i).folder, '\', str(i).name);
    output_file = strcat(str(i).folder, '\',str(i).name(1:end-6), '.nc') ;
    if ~exist(output_file,"file")
        command = sprintf('"%s" "%s" -netcdf "%s" >nul 2>&1', wgrib2_executable, file_path, output_file);      
        system(command);  %convert grib2 to netcdf:
    end

    % Get the nc file info 
    fileinfo = ncinfo(output_file);
    % get lon/lat data
    lon = ncread(output_file,'longitude'); % X
    lat = ncread(output_file,'latitude'); % Y
    lat = sort(lat,'descend')';

    % get time data
    time = fileinfo.Variables(3).Attributes(5).Value;
    datevec = datetime(time(1:end-3),"InputFormat","yyyy.MM.dd HH:mm:ss","TimeZone","UTC");
    datevec.Format = "dd-MMM-uuuu HH:mm:ss";
    S(i).latitude = lat;
    S(i).longitude = lon;
    S(i).time = datevec;
    S(i).significantWaveHeight = ncread(output_file,'HTSGW_surface');
    S(i).significantWavePeriod = ncread(output_file,'PERPW_surface');
    S(i).direction = ncread(output_file,'DIRPW_surface');
    S(i).windSpeed = ncread(output_file,'WIND_surface');
    S(i).windDirection = ncread(output_file,'WDIR_surface');
    delete(output_file)
end

struc = S;
end