data1 = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt\Sigma0_VH.hdr';
data2 = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt\Sigma0_VV.hdr';
%full_data = 'D:\UCT\EEE4022S\Data\CPT\Sigma0_VV.hdr';
smaller_subset_data = 'D:\UCT\EEE4022S\Data\CPT\largerSubset\Sigma0_VV.hdr';
fname = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt.h5';
str1 = readEnviHdr(data1);
str2 = readEnviHdr(data2);
%str3 = readEnviHdr(full_data);
str4 = readEnviHdr(smaller_subset_data);

% Number of Transects
n = 8;

data_VH = multibandread(str1.file,str1.size,str1.data_type,str1.header_offset,str1.interleave,str1.byte_order);
data_VH_norm = (data_VH - min(data_VH))./(max(data_VH) - min(data_VH));
data_VV = multibandread(str2.file,str2.size,str2.data_type,str2.header_offset,str2.interleave,str2.byte_order);
data_VV_norm = (data_VV - min(data_VV))./(max(data_VV) - min(data_VV));
%data_VV_full = multibandread(str3.file,str3.size,str3.data_type,str3.header_offset,str3.interleave,str3.byte_order);
data_VV_larger = multibandread(str4.file,str4.size,str4.data_type,str4.header_offset,str4.interleave,str4.byte_order);
h5_data = h5info(fname);

fileLoc = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt.h5';
testOut = getMetadataH5(fileLoc, 'abstracted');
testOutAtr = testOut.Attributes;

attribute_names = ["MISSION","SWATH", "BEAMS","PRODUCT", "orbit_cycle", "ABS_ORBIT"];
attribute_names2 = ["MISSION","SWATH", "BEAMS", "ABS_ORBIT"];
req_atributes = ["MISSION","orbit_cycle","first_line_time","antenna_pointing","PASS","centre_heading","slant_range_to_first_pixel","centre_lat","centre_lon","total_size"];
test_output = filterAttributesH5(testOutAtr,req_atributes);

test_meta_val = getAttributeValH5(test_output,req_atributes);
test_meta_val_for_inv = readMetadata(testOutAtr);
%transectData = zeros(512,512,2);
[transectData, startPos] = get512Transects(data_VV_larger,1,1,20,n);
[transectData2, startPos2] = get512Transects(data_VV_larger,1,1,45,n);
%transectData = data_VV(1:512,1:512);

%diference_1 = transectData2(:,:,2) - transectData(:,:,2);
%difference_2 = transectData2(:,:,3) - transectData(:,:,3);
test_ones = ones(1025,1025);
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
Ty=[0:size(tielat,2)-1]*stp(2)+1;
Tx=[0:size(tielat,1)-1]*stp(1)+1;
Iy=istart(2)+[0:size(subimg,2)-1]*strd(2);
Ix=istart(1)+[0:size(subimg,1)-1]*strd(1);
sublat=interp2(Ty',Tx,tielat,Iy',Ix);
sublon=interp2(Ty',Tx,tielon,Iy',Ix);

% Now make the map

m_proj('lambert','lon',[-34-45/60 -34-22/60],'lat',[16+34/60 16+44/60]);
m_pcolor(sublon,sublat,data_VV);shading flat;
m_grid('box','fancy','tickdir','out');
%m_ruler(1.03,[.15 .5],'ticklen',[.01]);
clim([0 2]);
colormap(gray);
title(titlestr)



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

figure(2)
subplot(2,3,1)
imshow(transectData(:,:,2));
title('VV Transect 2 20 deg')
subplot(2,3,2)
imshow(transectData(:,:,4));
title('VV Transect 4 20 deg')
subplot(2,3,3)
imshow(transectData(:,:,6));
title('VV Transect 6 20 deg')
subplot(2,3,4)
imshow(transectData2(:,:,2));
title('VV Transect 2 45 deg')
subplot(2,3,5)
imshow(transectData2(:,:,4));
title('VV Transect 4 45 deg')
subplot(2,3,6)
imshow(transectData2(:,:,6));
title('VV Transect 6 45 deg')
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
