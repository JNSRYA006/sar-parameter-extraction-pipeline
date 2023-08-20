data1 = 'D:\UCT\EEE4022S\Data\CPT\subset_1_cpt\Sigma0_VH.hdr';
data2 = 'D:\UCT\EEE4022S\Data\CPT\subset_1_cpt\Sigma0_VV.hdr';
str1 = readEnviHdr(data1);
str2 = readEnviHdr(data2);

data_VH = multibandread(str1.file,str1.size,str1.data_type,str1.header_offset,str1.interleave,str1.byte_order);
data_VH_norm = (data_VH - min(data_VH))./(max(data_VH) - min(data_VH));
data_VV = multibandread(str2.file,str2.size,str2.data_type,str2.header_offset,str2.interleave,str2.byte_order);
data_VV_norm = (data_VV - min(data_VV))./(max(data_VV) - min(data_VV));

%% Plots 
% figure(1)
% subplot(1,2,1)
% imshow(data_VH)
% title('VH Subset')
% subplot(1,2,2)
% imshow(data_VV)
% title('VV Subset')

figure(2)
subplot(2,2,1)
imshow(data_VH)
title('VH Subset')
subplot(2,2,2)
imshow(data_VH_norm)
title('VH Subset Normalised')
subplot(2,2,3)
imshow(data_VV)
title('VV Subset')
subplot(2,2,4)
imshow(data_VV_norm)
title('VV Subset Normalised')


% figure(3)
% subplot(1,2,1)
% imagesc(data_VH)
% title('VH Subset')
% subplot(1,2,2)
% imagesc(data_VV)
% title('VV Subset')
% %colormap jet
% %colorbar

figure(4)
subplot(1,2,1)
imagesc(data_VH)
title('VH Subset')
subplot(1,2,2)
imagesc(data_VH_norm)
title('VH Subset Normalised')
colormap jet
colorbar

figure(6)
subplot(1,2,1)
imagesc(data_VV)
title('VV Subset')
subplot(1,2,2)
imagesc(data_VV_norm)
title('VV Subset Normalised')
colormap jet
colorbar