function subWaveStruct = getSubsetWaveVals(waveStructFull,latBegin,lonBegin,latEnd,lonEnd)

    % Extract individual wave parameters from NetCDF file structure
    waveHeight = waveStructFull.significantWaveHeight;
    T = waveStructFull.significantWavePeriod;
    th = waveStructFull.direction; % degrees
    th = deg2rad(th); % rads
    lat = waveStructFull.latitude;
    lon = waveStructFull.longitude;
    lon = lon';
    time = waveStructFull.time;
    windSpeed = waveStructFull.windSpeed;
    windDirection = deg2rad(waveStructFull.windDirection);

    % Define and find the begining and end indices of the input lat and long values 
    latIndexBegin = find(round(lat,5) == latBegin);
    lonIndexBegin = find(round(lon,5) == lonBegin);
    latIndexEnd = find(round(lat,5) == latEnd);
    lonIndexEnd = find(round(lon,5) == lonEnd);

    % Extract the subset of lat and long values for all data types
    subWaveHeight = waveHeight(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd); % m
    subWavePeriod = T(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd); % s
    subWaveDirection = th(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd); % rads
    subLon = lon(lonIndexBegin:lonIndexEnd);
    subLat = lat(latIndexBegin:latIndexEnd);
    subWindSpeed = windSpeed(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd);
    subWindDirection = windDirection(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd);

    % Reconstruct structure with subset of data
    subWaveStruct.latitude = subLat;
    subWaveStruct.longitude = subLon;
    subWaveStruct.time = time;
    subWaveStruct.significantWaveHeight = subWaveHeight;
    subWaveStruct.significantWavePeriod = subWavePeriod;
    subWaveStruct.direction = subWaveDirection;
    subWaveStruct.windSpeed = subWindSpeed;
    subWaveStruct.windDirection = subWindDirection;

end