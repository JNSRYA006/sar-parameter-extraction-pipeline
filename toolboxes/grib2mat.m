function grib2mat(GribPath,SavePath)
% This function creats a structure file with the MRMS radar data and date
% time for each grib2 file in the input folder.

% -Inputs: - GribPath: folder path contain grib2 files
%          - SavePath: path to save the .mat structure
% -Output: - .mat file



str = dir(strcat(GribPath,'\*.grib2'));
S = struct('DirectionofWaves',[],'Latitude',[]);

% Define latitude and longitude ranges
lat_min = 33.97757;
lat_max = 34.22959;
lon_min = -106.9664;
lon_max = -106.6382;



for i = 1:length(str)
    % Define the directory where wgrib2 is located
    wgrib2_dir = 'C:\Users\ryanj\Downloads';
    
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

    % Extract values for the Pinos area:
    Values = ncread(output_file,fileinfo.Variables(1).Name);
    Values = flipud(Values);
    lat_mask = (lat >= lat_min) & (lat <= lat_max);
    lon_mask = (lon >= lon_min) & (lon <= lon_max);
    R = Values(lat_mask,lon_mask);

    time = fileinfo.Variables(3).Attributes(5).Value;
    datevec = datetime(time(1:end-3),"InputFormat","yyyy.MM.dd HH:mm:ss","TimeZone","UTC");
    datevec.Format = "dd-MMM-uuuu HH:mm:ss";
    S(i).RainData = fileinfo.Variables(1);
    S(i).Datetime = datevec;
    delete(output_file)
end
saveyear = num2str(S(1).Datetime.Year);
savemonth = num2str(S(1).Datetime.Month);
if length(savemonth)==1
    savemonth = strcat('0',savemonth);
end
saveday = num2str(S(1).Datetime.Day);
if length(savemonth)==1
    saveday = strcat('0',saveday);
end

savestr = strcat(SavePath,'\',saveyear,savemonth,saveday,'_MRMS');
save(savestr,'S')
end

