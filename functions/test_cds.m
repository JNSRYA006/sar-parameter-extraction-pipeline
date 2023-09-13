%% From Metadata Test
smaller_subset_data = 'D:\UCT\EEE4022S\Data\CPT\largerSubset\Sigma0_VV.hdr';
lat_long_subset = 'D:\UCT\EEE4022S\Data\CPT\lat_long_test\Sigma0_VV.hdr';
str2 = readEnviHdr(lat_long_subset);
str = readEnviHdr(smaller_subset_data);
%datasetName ="reanalysis-era5-single-levels";

% Number of Transects
n = 8;

data_VV = multibandread(str.file,str.size,str.data_type,str.header_offset,str.interleave,str.byte_order);
data_VV_2 = multibandread(str2.file,str2.size,str2.data_type,str2.header_offset,str2.interleave,str2.byte_order);

% Get metadata
fileLoc_H5 = 'D:\UCT\EEE4022S\Data\CPT\lat_long_test\lat_long_subset.h5';
testOut = getMetadataH5(fileLoc_H5, 'abstracted');
testOutAtr = testOut.Attributes;

attribute_names = ["first_line_time","first_near_lat", "first_near_long","first_far_lat", "first_far_long", "last_near_lat", "last_near_long","last_far_lat", "last_far_long"];
test_output = filterAttributesH5(testOutAtr,attribute_names);
test_meta_val = getAttributeValH5(test_output,attribute_names);
test_meta_val_for_inv = readMetadata(testOutAtr);
%% Break
% Get CDS Data from metadata
[U_10,V_10,long,lat,t,citation] = getWindVectorCDS(test_meta_val);

disp(citation);


U_10_t = permute(U_10, [2 1 3]);
V_10_t = permute(V_10, [2 1 3]);
%long_t = long.';
%lat_t = lat.';

% Plot data
[transectData, startPos] = get512Transects(data_VV_2,1,1,20,n);
[transectData2, startPos2] = get512Transects(data_VV_2,1,1,45,n);

%% World map plot
plotOnMapWind(long,lat,U_10_t,V_10_t,data_VV);
% worldmap('World')
% load coastlines
% plotm(coastlat,coastlon)


%% Plots
figure(1)
imshow(data_VV_2)
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
    annotate512Transect(startPos(i,1),startPos(i,2),i,'w','black',1);
    annotate512Transect(startPos2(i,1),startPos2(i,2),i,'r','r',0);
end

title('VV with transects and wind data from scene')

hold off

figure(2)
quiver(long,lat,U_10_t(:,:,1),V_10_t(:,:,1),'b','LineWidth',1)
hold on;
% Vector plot of second time frame
quiver(long,lat,U_10_t(:,:,2),V_10_t(:,:,2),'r','LineWidth',1)
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