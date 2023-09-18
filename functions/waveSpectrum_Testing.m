filePath = downloadNOAAWaveFile('https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20230905/18/wave/gridded/gfswave.t18z.gsouth.0p25.f300.grib2');
%% Real Deal
outStruct = getGribStruct("C:\Users\ryanj\Downloads","C:\Users\ryanj\OneDrive - University of Cape Town\4. Fourth Year\Second Semester\EEE4022S\repo\sar-parameter-extraction-pipeline\functions\");

%% Get variables
% Define constants
g = 9.81;

waveHeight = outStruct.significantWaveHeight;
T = outStruct.significantWavePeriod;
th = outStruct.direction; % degrees
th = deg2rad(th); % rads
lat = outStruct.latitude;
lon = outStruct.longitude.';
time = outStruct.time;

latBegin = -35;
lonBegin = 30;
latEnd = -40;
lonEnd = 35;

latIndexBegin = find(round(lat,5) == latBegin);
lonIndexBegin = find(round(lon,5) == lonBegin);
latIndexEnd = find(round(lat,5) == latEnd);
lonIndexEnd = find(round(lon,5) == lonEnd);


testWaveHeight = waveHeight(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd); % m
testWavePeriod = T(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd); % s
testWaveDirection = th(lonIndexBegin:lonIndexEnd,latIndexBegin:latIndexEnd); % rads
testLon = lon(lonIndexBegin:lonIndexEnd);
testLat = lat(latIndexBegin:latIndexEnd);

%% Calculate Wave Spectrum
lambda = g.*testWavePeriod.^2/2*pi;
k = (4*pi^2)./(g.*testWavePeriod);
k_x = k.*cos(testWaveDirection);
k_y = k.*sin(testWaveDirection);

% Calculate the FFT of wave height data
height_fft = fft2(testWaveHeight);

% Shift the zero frequency components to the center of the spectrum
height_fft_shifted = fftshift(height_fft);

% Calculate the magnitude (amplitude) spectrum
amplitude_spectrum = abs(height_fft_shifted);

% Create a grid of k_x and k_y values corresponding to the FFT output
[Ny, Nx] = size(testWaveHeight);
%k_x = fftshift(1:Nx) - floor(Nx/2);
%k_y = fftshift(1:Ny) - floor(Ny/2);

% Assuming you have sorted k_x and k_y
k_x_sorted = sort(k_x);
k_y_sorted = sort(k_y);

% Define grid spacing (adjust as needed)
dk_x = k_x_sorted(2) - k_x_sorted(1);
dk_y = k_y_sorted(2) - k_y_sorted(1);

% Plot the amplitude spectrum
figure;
contourf(k_x_sorted, k_y_sorted, amplitude_spectrum, 'LineStyle', 'none');
xlabel('k_x');
ylabel('k_y');
title('Amplitude Spectrum');
colorbar;

%% Plot
figure(1)

m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
'lon',[min(lon(:)) max(lon(:))])

% Next, plot the field using the M_MAP version of pcolor.

m_pcolor(testLon,testLat,k.');
shading flat;

% Add a coastline and axis values.

m_coast('patch',[.7 .7 .7])
m_grid('box','fancy')

% Add a colorbar and title.

colorbar
title('WAVEWATCH III Wavenumber from NOMADS');
