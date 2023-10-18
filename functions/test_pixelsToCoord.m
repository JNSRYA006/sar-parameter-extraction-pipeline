%% Test pixel values to be mapped to corresponding latitude and longitude
% Import data
filepath = "D:\UCT\EEE4022S\Data\CPT\smaller_subset_incidence.nc";
% Number of Transects
n = 2;

% Import data values
ncImport = ncinfo(filepath);
VV_nc = ncread(filepath,'Sigma0_VV');
[transectData_nc, startPos_nc] = get512Transects(VV_nc,1,1,20,n);
[transectData2_nc, startPos2] = get512Transects(VV_nc,1,1,45,n);

% Required Metadata
meta_nc = ncinfo(filepath,'metadata');
req_atributes = ["first_near_lat","first_near_long","first_far_lat","first_far_long","last_near_lat","last_near_long","last_far_lat","last_far_long","centre_lat","centre_lon", "num_output_lines","num_samples_per_line"];
meta_nc = filterAttributesNetCDF(meta_nc.Attributes, req_atributes);

%% Define start and end lat and long
% Want a 5 x 2 matrix of lat and long
% Row 1 - first_near
% Row 2 - first_far
% Row 3 - last_near
% Row 4 - last_far
% Row 5 - Centre
coords = zeros(5,2);
coords(1,:) = [meta_nc(1).Value,meta_nc(2).Value];
coords(2,:) = [meta_nc(3).Value,meta_nc(4).Value];
coords(3,:) = [meta_nc(5).Value,meta_nc(6).Value];
coords(4,:) = [meta_nc(7).Value,meta_nc(8).Value];
coords(5,:) = [meta_nc(9).Value,meta_nc(10).Value];

%coords_order = sort(coords(1:4,:), 1, 'descend');
max_lat = max(coords(1:4,1));
max_lon = max(coords(1:4,2));
min_lat = min(coords(1:4,1));
min_lon = min(coords(1:4,2));
%% Define size 
% 1 x 2 matrix
size = zeros(1,2);
size(1,:) = [meta_nc(12).Value,meta_nc(11).Value];

%% Create linspace of lat and long over number of pixels
lat_interpolated = linspace(max_lat,min_lat,size(1,1));
lon_interpolated = linspace(min_lon,max_lon,size(1,2));

%% Mmap plot
figure(3)
m_proj('lambert','long', [min_lon max_lon],'lat',[min_lat max_lat]);
clf;
m_pcolor(lon_interpolated,lat_interpolated,VV_nc); shading flat;
%imshow(VV_nc); set(gca, 'YDir','normal');
%m_vec(1,long,lat,U_10_t(:,:,1),V_10_t(:,:,1),'k','shaftwidth',0.8);
m_grid('box','fancy','tickdir','out');
m_ruler(1.03,[.15 .5],'ticklen',.01);
%colormap(gray);
%m_grid('tickdir','out','linewi',3,'fontsize',14);
title('Test Plot using mmap');


%% Plot SAR data on lat lon grid instead of pixels grid

lonOffset = min_lon;
latOffset = min_lat;
lonScale = (max_lon - min_lon)/size(1,2);
latScale = (max_lat - min_lat)/size(1,1);

figure(1)
imshow(VV_nc)

grid on;
hold on;

for i = 1:n
    annotate512Transect(startPos_nc(i,1),startPos_nc(i,2),i,'w','black',1);
    annotate512Transect(startPos2(i,1),startPos2(i,2),i,'r','r',0);
end
xMin = (min_lon - lonOffset) / lonScale; % Convert lonMin to pixel x-coordinate
xMax = (max_lon - lonOffset) / lonScale; % Convert lonMax to pixel x-coordinate
yMin = (min_lat - latOffset) / latScale; % Convert latMin to pixel y-coordinate
yMax = (max_lat - latOffset) / latScale; % Convert latMax to pixel y-coordinate

title('VV with transects')
% Define custom axes for the image
xlim([xMin, xMax]);
ylim([yMin, yMax]);

% Manually add axis labels for longitude and latitude
xTicks = linspace(1, size(1,2), 5);  % Adjust as needed
xLabels = linspace(min_lon, max_lon, 5);        % Adjust as needed
yTicks = linspace(1, size(1,1), 5); % Adjust as needed
yLabels = linspace(min_lat, max_lat, 5);        % Adjust as needed
xticks(xTicks);
yticks(yTicks);
% set(gca, 'XTick', xTicks, 'XTickLabel', xLabels);
% set(gca, 'YTick', yTicks, 'YTickLabel', yLabels);

% Optionally, label the axes
xlabel('Longitude');
ylabel('Latitude');
hold off

%% Subplots of transects
figure(2)
subplot(1,2,1)
imshow(transectData_nc(:,:,1))
title('VH Transect 1')
subplot(1,2,2)
imshow(transectData_nc(:,:,2))
title('VH Transect 2')
% subplot(3,2,3)
% imshow(transectData_nc(:,:,3))
% title('VH Transect 3')
% subplot(3,2,4)
% imshow(transectData_nc(:,:,4))
% title('VH Transect 4')
% subplot(3,2,5)
% imshow(transectData_nc(:,:,5))
% title('VV Transect 5')
% subplot(3,2,6)
% imshow(transectData_nc(:,:,6))
% title('VV Transect 6')