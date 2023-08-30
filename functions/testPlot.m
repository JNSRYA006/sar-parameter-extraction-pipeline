data1 = 'D:\UCT\EEE4022S\Data\CPT\subset_1_cpt\Sigma0_VH.hdr';
data2 = 'D:\UCT\EEE4022S\Data\CPT\subset_1_cpt\Sigma0_VV.hdr';
%full_data = 'D:\UCT\EEE4022S\Data\CPT\Sigma0_VH.hdr';
str1 = readEnviHdr(data1);
str2 = readEnviHdr(data2);
%str3 = readEnviHdr(full_data);

% Number of Transects
n = 2;

data_VH = multibandread(str1.file,str1.size,str1.data_type,str1.header_offset,str1.interleave,str1.byte_order);
data_VH_norm = (data_VH - min(data_VH))./(max(data_VH) - min(data_VH));
data_VV = multibandread(str2.file,str2.size,str2.data_type,str2.header_offset,str2.interleave,str2.byte_order);
data_VV_norm = (data_VV - min(data_VV))./(max(data_VV) - min(data_VV));
%data_VH_full = multibandread(str3.file,str3.size,str3.data_type,str3.header_offset,str3.interleave,str3.byte_order);

fileLoc = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt.h5';
testOut = getMetadataH5(fileLoc, 'abstracted');
testOutAtr = testOut.Attributes;

attribute_names = ["MISSION","SWATH", "BEAMS","PRODUCT", "orbit_cycle", "ABS_ORBIT"];
attribute_names2 = ["MISSION","SWATH", "BEAMS", "ABS_ORBIT"];
test_output = filterAttributesH5(testOutAtr,attribute_names);

test_meta_val = getAttributeValH5(test_output,attribute_names2);
%transectData = zeros(512,512,2);
[transectData, startPos] = get512Transects(data_VV,500,500,45,n);
[transectData2, startPos2] = get512Transects(data_VH_full,1,1,30,n);
%transectData = data_VV(1:512,1:512);

% diference_1 = transectData2(:,:,2) - transectData(:,:,2);
% difference_2 = transectData2(:,:,3) - transectData(:,:,3);

%% Plots 

figure(1)
% subplot(3,2,1)
imshow(data_VV)
hold on;

for i = 1:n
    annotate512Transect(startPos(i,1),startPos(i,2),i,'w','black',1);
    annotate512Transect(startPos2(i,1),startPos2(i,2),i,'r','r',0);
end

title('VV with transects')
hold off

figure(2)
subplot(2,3,1)
imshow(transectData(:,:,1));
title('VV Transect 1 45 deg')
subplot(2,3,2)
imshow(transectData(:,:,2));
title('VV Transect 2 45 deg')
subplot(2,3,3)
imshow(transectData(:,:,3));
title('VV Transect 3 45 deg')
subplot(2,3,4)
imshow(transectData(:,:,1));
title('VV Transect 1 30 deg')
subplot(2,3,5)
imshow(transectData(:,:,2));
title('VV Transect 2 30 deg')
subplot(2,3,6)
imshow(transectData(:,:,3));
title('VV Transect 3 30 deg')
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
% imshow(diference_1)
% title('Difference 1')
% subplot(2,2,4)
% imshow(difference_2)
% title('Difference 2')

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
% imagesc(data_VH)
% title('VH Subset')
% subplot(1,2,2)
% imagesc(data_VV)
% title('VV Subset')
% %colormap jet
% %colorbar

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
