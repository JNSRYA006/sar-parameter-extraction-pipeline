grib='akq_nwps_CG0_Trkng_20230909_1200.grib2';

nco=ncgeodataset(grib);
nco.variables
% ans = 
%     'Wind_direction_from_which_blowing_degree_true_surface'
%     'Wind_speed_surface'
%     'u-component_of_wind_surface'
%     'v-component_of_wind_surface'
%     'Significant_height_of_combined_wind_waves_and_swell_surface'
%     'Direction_of_wind_waves_degree_true_surface'
%     'Significant_height_of_wind_waves_surface'
%     'Mean_period_of_wind_waves_surface'
%     'Direction_of_swell_waves_degree_true_ordered_sequence_of_data'
%     'Significant_height_of_swell_waves_ordered_sequence_of_data'
%     'Mean_period_of_swell_waves_ordered_sequence_of_data'
%     'Primary_wave_direction_degree_true_surface'
%     'Primary_wave_mean_period_surface'
%     'lat'
%     'lon'
%     'ordered_sequence_of_data'
%     'time'

% Extract the Significant Height of Swell Waves field.

param='Significant_height_of_swell_waves_ordered_sequence_of_data';
waveheight=nco{param}(1,1,:,:);
lat=nco{'lat'}(:);
lon=nco{'lon'}(:);

% From this point on the code is identical to the previous example:

waveheight=double(squeeze(waveheight));
lat=double(lat);
lon=double(lon);
m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
  'lon',[min(lon(:)) max(lon(:))])
m_pcolor(lon,lat,waveheight);
shading flat;
m_coast('patch',[.7 .7 .7]);
m_grid('box','fancy')
colorbar
title('Example 2: WAVEWATCH III Significant Wave Height from GRiB');