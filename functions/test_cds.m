%% From Metadata Test
% Import data
filepath = "D:\UCT\EEE4022S\Data\CPT\small_subset.nc";
% Number of Transects
n = 2;

% Import data values
ncImport = ncinfo(filepath);
VV_nc = ncread(filepath,'Sigma0_VV');
[transectData_nc, startPos_nc] = get512Transects(VV_nc,1,1,20,n);
[transectData2_nc, startPos2] = get512Transects(VV_nc,1,1,45,n);

% Required Metadata
meta_nc = ncinfo(filepath,'metadata');
req_atributes = ["first_near_lat","first_near_long","first_far_lat","first_far_long","last_near_lat","last_near_long","last_far_lat","last_far_long","centre_lat","centre_lon", "num_output_lines","num_samples_per_line","first_line_time"];
meta_nc = filterAttributesNetCDF(meta_nc.Attributes, req_atributes);
%% Break
% Get CDS Data from metadata
[U_10,V_10,long,lat,t,citation,directional,peakedness] = getWindVectorCDS(filepath);

disp(citation);


U_10_t = permute(U_10, [2 1 3]);
V_10_t = permute(V_10, [2 1 3]);
%long_t = long.';
%lat_t = lat.';


%% Test Heading plot
% Single heading value (in degrees)
heading = deg2rad(90 - 344.3488); % Replace with your desired heading angle

% Create a unit vector based on the heading angle
u = cos(heading); % Calculate x-component of the unit vector
v = sin(heading); % Calculate y-component of the unit vector

%% Mmap plot
% figure(1)
% 
% m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
% 'lon',[min(lon(:)) max(lon(:))])
% 
% Next, plot the field using the M_MAP version of pcolor.
% 
% m_pcolor(lon,lat,waveHeight.');
% shading flat;
% 
% Add a coastline and axis values.
% 
% m_coast('patch',[.7 .7 .7])
% m_grid('box','fancy')
% 
% ax.YLabel.String='Wind speeds m/s';
% 
% 
% Add a colorbar and title.
% 
% colorbar
% title('WAVEWATCH III Wavenumber from NOMADS');

%% Plots
figure(1)
imshow(VV_nc)
h = gca;
h.Visible = 'On';
hold on;

% Vector plot of first time frame
wind_time_1 = quiver(long,lat,U_10_t(:,:,1),V_10_t(:,:,1),'b','LineWidth',1);
% Vector plot of second time frame
wind_time_2 = quiver(long,lat,U_10_t(:,:,2),V_10_t(:,:,2),'r','LineWidth',1);
legend([wind_time_1, wind_time_2], ['Wind at: ', '17:00'], ['Wind at: ', '18:00']);
% Plot rectangles of transects
for i = 1:n
    annotate512Transect(startPos_nc(i,1),startPos_nc(i,2),i,'w','black',1);
    annotate512Transect(startPos2(i,1),startPos2(i,2),i,'r','r',0);
end

title('VV with transects and wind data from scene')

hold off

figure(2)
quiver(long,lat,U_10_t(:,:,1),V_10_t(:,:,1),'b','LineWidth',1)
hold on;
% Vector plot of second time frame
quiver(long,lat,U_10_t(:,:,2),V_10_t(:,:,2),'r','LineWidth',1)
%quiver(17.6424, -34.0603, u, v, 'bl', 'LineWidth', 2); % Assuming origin at (centre_lon, centre_lat)
hold off;

%% Hard coded test
datasetName ="reanalysis-era5-single-levels";

% Choose variables
options.product_type = "reanalysis";
options.format = "netcdf";
options.grid = ["1.0","1.0"];
%options.version = "1_0";
options.variable = ["10m_u_component_of_wind","10m_v_component_of_wind"];
options.year = "2023"; 
options.month = "09";
options.day = "06";
options.time = ["17:00","01:00"];
%options.time = "17:00";
options.area = ["-32","14","-36","20"];



% Download data
[downloadedFilePaths,citation] = climateDataStoreDownload(datasetName,options);

% Display citation
disp(citation);


U_10 = ncread(downloadedFilePaths,'u10');
V_10 = ncread(downloadedFilePaths,'v10');
long = ncread(downloadedFilePaths,'longitude');
lat = ncread(downloadedFilePaths,'latitude');
t = ncread(downloadedFilePaths,'time');

U_10_t = permute(U_10, [2 1 3]);
V_10_t = permute(V_10, [2 1 3]);
long_t = long.';
lat_t = lat.';

figure(3)
% Vector plot of first time frame
quiver(long,lat,U_10_t(:,:,1),V_10_t(:,:,1),'b','LineWidth',1)
hold on;
% Vector plot of second time frame
quiver(long,lat,U_10_t(:,:,2),V_10_t(:,:,2),'r','LineWidth',1)
hold off;

%Clean up the downloaded files
% delete(downloadedFilePaths)