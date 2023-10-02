function [U_10,V_10,long,lat,t,citation,directionalWidth,peakedness] = getWindVectorCDS(metadataFilePath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

attribute_names = ["first_line_time","first_near_lat", "first_near_long","first_far_lat", "first_far_long", "last_near_lat", "last_near_long","last_far_lat", "last_far_long"];
datasetName ="reanalysis-era5-single-levels";

% Get CDS Data from metadata
% Format data appropriately
meta_nc = ncinfo(metadataFilePath,'metadata');
meta_nc = filterAttributesNetCDF(meta_nc.Attributes, attribute_names);

area_meta = formatLatLonValues(meta_nc);

% Get datetime data
% datetime_meta = datetime(metadata(strcmp({metadata.Name},attribute_names(1))).Value, 'Format', 'yyyy-MM-dd HH:mm:ss');
% hour_meta = string([num2str(hour(datetime_meta)), ':00']);
% next_hour_meta = string([num2str(hour(datetime_meta) + 1), ':00']);
% time_meta = [hour_meta, next_hour_meta];
% day_meta = string(num2str(day(datetime_meta)));
% month_meta = string(num2str(month(datetime_meta)));
% year_meta = string(num2str(year(datetime_meta)));

[time_meta, day_meta, month_meta, year_meta] = formatDateValues(meta_nc,attribute_names);

% Set option for api call
options.product_type = "reanalysis";
options.format = "netcdf";
options.grid = ["1.0","1.0"];
options.variable = ["10m_u_component_of_wind","10m_v_component_of_wind", "wave_spectral_directional_width_for_wind_waves", "wave_spectral_peakedness"];
options.year = year_meta; 
options.month = month_meta;
options.day = day_meta;
options.time = time_meta;
options.area = area_meta;



% Download netCDF data and citation
[downloadedFilePaths,citation] = climateDataStoreDownload(datasetName,options);
fileStruct = ncinfo(downloadedFilePaths);
% Extract required parameters
U_10 = ncread(downloadedFilePaths,'u10');
V_10 = ncread(downloadedFilePaths,'v10');
long = ncread(downloadedFilePaths,'longitude');
lat = ncread(downloadedFilePaths,'latitude');
t = ncread(downloadedFilePaths,'time');
directionalWidth = ncread(downloadedFilePaths,'dwww');
peakedness = ncread(downloadedFilePaths,'wsp');

%Clean up the downloaded files
delete(downloadedFilePaths)

end

function area = formatLatLonValues(metadata)

% Get all lat values
lat_near1 = metadata(2).Value;
lat_far1 = metadata(4).Value;
lat_near2 = metadata(6).Value;
lat_far2 = metadata(8).Value;
lat = [lat_near1,lat_near2,lat_far1,lat_far2];

% Get all long values
lon_near1 = metadata(3).Value;
lon_far1 = metadata(5).Value;
lon_near2 = metadata(7).Value;
lon_far2 = metadata(9).Value;
lon = [lon_near1,lon_near2,lon_far1,lon_far2];

% Determine max and min values of latitude and longitude and set coordinate 
% values to whole numbers to keep CDS api call happy
lat_N = string(num2str(ceil(max(lat))));
lat_S = string(num2str(floor(min(lat))));
lon_E = string(num2str(ceil(min(lon))));
lon_W = string(num2str(floor(max(lon))));

area = [lat_N,lon_E,lat_S,lon_W];

end

function [times, day_str, month_str, year_str] = formatDateValues(metadata,attribute_names)

% Get datetime data
datetime_meta = datetime(metadata(1).Value, 'Format', 'yyyy-MM-dd HH:mm:ss');
hour_meta = string([num2str(hour(datetime_meta)), ':00']);
next_hour_meta = string([num2str(hour(datetime_meta) + 1), ':00']);
times = [hour_meta, next_hour_meta];
day_str = string(num2str(day(datetime_meta)));
month_str = string(num2str(month(datetime_meta)));
year_str = string(num2str(year(datetime_meta)));

end
