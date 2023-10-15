%grib2mat("C:\Users\ryanj\OneDrive - University of Cape Town\4. Fourth Year\Second Semester\EEE4022S\repo\sar-parameter-extraction-pipeline\toolboxes\","C:\Users\ryanj\Downloads");

outStruct = getGribStruct("C:\Users\ryanj\Downloads","C:\Users\ryanj\OneDrive - University of Cape Town\4. Fourth Year\Second Semester\EEE4022S\repo\sar-parameter-extraction-pipeline\toolboxes\");
%%
filePath = downloadNOAAWaveFile('https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.20230913/12/wave/gridded/gdaswave.t12z.gsouth.0p25.f000.grib2');
%% Real Deal
outStruct = getGribStruct("C:\Users\ryanj\Downloads","C:\Users\ryanj\OneDrive - University of Cape Town\4. Fourth Year\Second Semester\EEE4022S\repo\sar-parameter-extraction-pipeline\functions\");
%%
waveheight = outStruct.significantWaveHeight;
wavePeriod = outStruct.significantWavePeriod;
waveDirection = deg2rad(outStruct.direction);
lat = outStruct.latitude;
lon = outStruct.longitude.';
heading = 334.3488;

waveheight=double(squeeze(waveheight)).';
wavePeriod=double(squeeze(wavePeriod)).';
waveDirection=double(squeeze(waveDirection)).';

% Plot the field using M_MAP.  Start with setting the map
% projection using the limits of the lat/lon data itself:
figure(1)
subplot(3,1,1)
m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
'lon',[min(lon(:)) max(lon(:))])

% Next, plot the field using the M_MAP version of pcolor.

m_pcolor(lon,lat,waveheight);
shading flat;

% Add a coastline and axis values.

m_coast('patch',[.7 .7 .7])
m_grid('box','fancy')

% Add a colorbar and title.

colorbar
title('WAVEWATCH III Significant Wave Height from NOMADS');


% Plot the field using M_MAP.  Start with setting the map
% projection using the limits of the lat/lon data itself:
subplot(3,1,2)
m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
'lon',[min(lon(:)) max(lon(:))])

% Next, plot the field using the M_MAP version of pcolor.

m_pcolor(lon,lat,wavePeriod);
shading flat;

% Add a coastline and axis values.

m_coast('patch',[.7 .7 .7])
m_grid('box','fancy')

% Add a colorbar and title.

colorbar
title('WAVEWATCH III Significant Wave Period from NOMADS');

subplot(3,1,3)
m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
'lon',[min(lon(:)) max(lon(:))])

% Next, plot the field using the M_MAP version of pcolor.

m_pcolor(lon,lat,waveDirection);
shading flat;

% Add a coastline and axis values.

m_coast('patch',[.7 .7 .7])
m_grid('box','fancy')

% Add a colorbar and title.

colorbar
title('WAVEWATCH III Wave Direction from NOMADS');