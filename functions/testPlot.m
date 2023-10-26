%% nc import data test
% Number of Transects
n = 3;

% Test getting values
ncImport = ncinfo('D:\UCT\EEE4022S\Data\CPT\landSeaCP.nc');
VV_nc = ncread('D:\UCT\EEE4022S\Data\CPT\landSeaCP.nc','Sigma0_VV');
%%
[transectData_nc, startPos_nc] = get512Transects(VV_nc,1,1,30,n);
[transectData_th, startPos_th] = get512Transects(th,1,1,30,n);
%[transectData2, startPos2] = get512Transects(data_VV_larger,1,1,45,n);

%% Test metadata
meta_nc = ncinfo('D:\UCT\EEE4022S\Data\CPT\larger_subset.nc','metadata');
req_atributes = ["MISSION","orbit_cycle","first_line_time","antenna_pointing","PASS","centre_heading","slant_range_to_first_pixel","centre_lat","centre_lon","total_size"];
%req_atributes = ["MISSION","SWATH", "BEAMS", "ABS_ORBIT"];
meta_nc = filterAttributesNetCDF(meta_nc.Attributes, req_atributes);
%% RMSE
transectRMSE_1_2 = sum(rmse(transectData_nc(:,:,1),transectData_nc(:,:,2)));
transectRMSE_1_3 = sum(rmse(transectData_nc(:,:,1),transectData_nc(:,:,3)));
transectRMSE_2_3 = sum(rmse(transectData_nc(:,:,2),transectData_nc(:,:,3)));
transectRMSE_same = sum(rmse(transectData_nc(:,:,1),transectData_nc(:,:,1))); 
transectMSE_1_2 = transectRMSE_1_2^2;
transectMSE_1_3 = transectRMSE_1_3^2;
transectMSE_2_3 = transectRMSE_2_3^2;
transectMSE_1_1 = transectRMSE_same^2;

transectMSEStruct = transectMSE(transectData_nc);
%% GPT Sar Spectra
% Calculate the 2D Fourier Transform
spectrum = fftshift(fft2(transectData));
spectrum = abs(spectrum);
spectrum_index = spectrum > 20;
spectrum = spectrum_index.*spectrum;

% Calculate the spatial frequencies
[M, N] = size(transectData(:,:,1));
dx = spectral_bw; % Change this value to match your data scale
dy = spectral_bw; % Change this value to match your data scale
% kx = (-M/2:M/2-1) / (M * dx);
% ky = (-N/2:N/2-1) / (N * dy);
kx = k_x;
ky = k_y;

% Calculate Δk (the spacing between k values)
dkx = kx(2) - kx(1);
dky = ky(2) - ky(1);

% Define the range of indices to crop spectrum
startIndexX = round(M/2) - 49; % Adjust as needed
endIndexX = round(M/2) + 50; % Adjust as needed
startIndexY = round(N/2) - 49; % Adjust as needed
endIndexY = round(N/2) + 50; % Adjust as needed

% Select the subset of data for the middle 100x100 points
spectrumSubset = spectrum(startIndexX:endIndexX, startIndexY:endIndexY);
kxSubset = kx(startIndexX:endIndexX);
kySubset = ky(startIndexY:endIndexY);

% Create a grid for the contour plot
%[X, Y] = meshgrid(kxSubset, kySubset);
[X, Y] = meshgrid(kx, ky);

% Plot the SAR spectrum with k_x/Δk and k_y/Δk axes using contourf
figure;
%contourf(X, Y, 20*log10(abs(spectrum)), 20, 'LineColor', 'none');
contour(X, Y, 20*log10(spectrum));
colormap('jet');
colorbar;
yline(0);
xline(0);

% % Set the X and Y axis ticks and labels
% xticks(unique([min(kxSubset):1:max(kxSubset)]));
% yticks(unique([min(kySubset):1:max(kySubset)]));
% xticklabels(unique(string(kxSubset)));
% yticklabels(unique(string(kySubset)));

xlabel('k_x/Δk');
ylabel('k_y/Δk');

title('SAR Spectrum with k_x/Δk and k_y/Δk Axes');


%% Transect SAR spectrum with circle plots
%% Get SAR intensity power spectrum
% Calculate the 2D Fourier transform of the intensity data
intensityFFT = abs(fftshift(fft2(transectData)));
dk_x = 0.005890997229533;
dk_y = -0.005226593472072;
%%
%contourf(k_x./dk_x,k_y./dk_y,20*log10(intensityFFT(:,:,1)))
%contourf(20*log10(intensityFFT(:,:,1)))
%surf(20*log10(intensityFFT),'LineStyle','none')
%imagesc(20*log10(intensityFFT))
figure(3);
contourf(20*log10(intensityFFT(:,:,1)))
colorbar;
title(['Transect 1 SAR Spectrum at ',num2str(th),' degrees'])
hold on;
% Define the center and radius of the circle
centerX = 512/2;  % X-coordinate of the center
centerY = 512/2;  % Y-coordinate of the center
radii = [50,100,200,400];  % Radius of the circle
radii = abs(fft(radii));
radiiCol = ['r','w','b','bl'];
for i=1:length(radii)
    radius = radii(i);
    rectangle('Position', [centerX - radius,centerY - radius, 2*radius, 2*radius], ...
    'Curvature', [1, 1], 'EdgeColor', radiiCol(i));  % 'Curvature' set to [1, 1] makes it a circle
end
axis equal;  % Equal aspect ratio for proper circle visualization
hold off;
%%
figure(4);
contourf(20*log10(intensityFFT(:,:,2)))
colorbar;
title(['Transect 2 SAR Spectrum at ',num2str(th),' degrees'])
figure(5);
contourf(20*log10(intensityFFT(:,:,3)))
colorbar;
title(['Transect 3 SAR Spectrum at ',num2str(th),' degrees'])







%%
figure;
imshow(VV_nc)
%title('S1A VV data from 28-Jul-2023 at 34.78S,16.77E with transects at 30$\degree$ shown',Interpreter='latex');
hold on;

for i = 1:n
    annotate512Transect(startPos_nc(i,1),startPos_nc(i,2),i,'red','black',1);
end

hold off
% 
figure;
subplot(1,3,1)
imshow(transectData_nc(:,:,1));
%title('VV Transect 1 30 deg')
%figure;
subplot(1,3,2)
imshow(transectData_nc(:,:,2));
%title('VV Transect 2 30 deg')
%figure;
subplot(1,3,3)
imshow(transectData_nc(:,:,3));
%title('VV Transect 3 30 deg')

%% Plot - M-Map (SAR)
titlestr='Test SAR Plot with m-maps';
datsize=double([1025 1025]);

tielat=h5read(fname,'/tie_point_grids/latitude');
tielon=h5read(fname,'/tie_point_grids/longitude');
stp=[h5readatt(fname,'/tie_point_grids/latitude','sub_sampling_x') ...
     h5readatt(fname,'/tie_point_grids/latitude','sub_sampling_y') ];

subf = filter2(ones(3,3)/9,data_VV);

% Now generate lat/lon for all pixels by interpolating from
% the tie points.
% Ty=[0:size(tielat,2)-1]*stp(2)+1;
% Tx=[0:size(tielat,1)-1]*stp(1)+1;
% Iy=istart(2)+[0:size(subimg,2)-1]*strd(2);
% Ix=istart(1)+[0:size(subimg,1)-1]*strd(1);
% sublat=interp2(Ty',Tx,tielat,Iy',Ix);
% sublon=interp2(Ty',Tx,tielon,Iy',Ix);

% Now make the map
% 
% m_proj('lambert','lon',[-34-45/60 -34-22/60],'lat',[16+34/60 16+44/60]);
% m_pcolor(sublon,sublat,data_VV);shading flat;
% m_grid('box','fancy','tickdir','out');
% %m_ruler(1.03,[.15 .5],'ticklen',[.01]);
% clim([0 2]);
% colormap(gray);
% title(titlestr)



%% Plots - M_map
m_proj('ortho','lat',-34.6305','long',16.6026);
%m_plot(16.6026,-34.6305,data_VH);
m_grid('linestyle','-','xticklabels',[],'yticklabels',[],'ytick',[-34:16:80]);
xlabel('Orthographic Projection','visible','on');


%% Plots - MATLAB

% figure(1)
% % subplot(3,2,1)
% imshow(data_VV_larger)
% hold on;
% 
% for i = 1:n
%     annotate512Transect(startPos(i,1),startPos(i,2),i,'w','black',1);
%     annotate512Transect(startPos2(i,1),startPos2(i,2),i,'r','r',0);
% end
% 
% title('VV with transects')
% hold off

% subplot(3,2,2)
% imshow(transectData(:,:,1))
% title('VH Transect 1')
% subplot(3,2,3)
% imshow(transectData(:,:,2))
% title('VH Transect 2')
% subplot(3,2,4)
% imshow(transectData(:,:,3))
% title('VH Transect 3')
% subplot(3,2,5)
% imshow(transectData(:,:,4))
% title('VH Transect 4')
% subplot(3,2,6)
% imshow(transectData2(:,:,5))
% title('VV Transect 5')
% subplot(2,2,1)
% imshow(transectData(:,:,2))
% title('VV Transect 2')
% subplot(2,2,2)
% imshow(transectData(:,:,3))
% title('VV Transect 3')
% subplot(2,2,3)

figure(3)
subplot(1,2,1)
imshow(data_VH)
title('VH Data')
subplot(1,2,2)
imshow(data_VV)
title('VV Data')

figure(4)
pcolor(data_VV);
colorbar;
shading interp
lims = clim;
max(data_VV);
% figure(1)
% subplot(1,2,1)
% imagesc(data_VH)
% title('VH Subset')
% subplot(1,2,2)
% imagesc(data_VV)
% title('VV Subset')

% figure(2)
% subplot(2,2,1)
% imshow(data_VH)
% hold on;
% axis on;
% text(200,256,'1',FontSize=10,Color='w');
% rectangle('Position',[1,1,512,512],'EdgeColor', 'w','LineWidth', 1,'LineStyle','-')
% title('VH Subset')
% hold off
% subplot(2,2,2)
% imshow(data_VH_norm)
% title('VH Subset Normalised')
% subplot(2,2,3)
% imshow(data_VV)
% title('VV Subset')
% subplot(2,2,4)
% imshow(data_VV_norm)
% title('VV Subset Normalised')


% figure(3)
% subplot(1,2,1)
% imagesc(transectData(:,:,1))
% title('VH Subset')
% subplot(1,2,2)
% imagesc(transectData(:,:,2))
% title('VV Subset')
% colormap jet
% colorbar

% figure(4)
% subplot(1,2,1)
% imagesc(data_VH)
% title('VH Subset')
% subplot(1,2,2)
% imagesc(data_VH_norm)
% title('VH Subset Normalised')
% colormap jet
% colorbar

% figure(6)
% subplot(1,2,1)
% imagesc(data_VV)
% title('VV Subset')
% subplot(1,2,2)
% imagesc(data_VV_norm)
% title('VV Subset Normalised')
% colormap jet
% colorbar
