%% Required files
% Import data
filepath = "D:\UCT\EEE4022S\Data\CPT\smaller_subset_incidence.nc";
%[Filename,filepath]=uigetfile('*.nc','Choose the exported .nc file from SNAP');
%vfilen=[filepath,Filename];
% Import data values
ncImport = ncinfo(filepath);
VV_nc = ncread(filepath,'Sigma0_VV');

%% Required Metadata
metadata = ncinfo(filepath,'metadata');
th = ncread(filepath,'Incidence_Angle');
r = ones(2001);
%% Testing orbital velocity covariance function
[f_v_r] = orbitalVelocityCovariance(k,k_y,E_k,metadata,th,r);
%% Testing autocovariance function of RAR image intensity
[f_r_r] = rarImageIntensityAutocovariance(k,k_y,E_k,metadata,th,r);
%% Testing covariance function of RAR image intensity and NL velocity
[f_rv_r] = rarImageIntensityCovariance(k,k_y,E_k,metadata,th,r);
%% Quasi-linear Approximation
[p_s_coeff,beta] = quasilinearCoeff(k,k_y,k_x,E_k,metadata,th);
%% Calculate Power Spectrum Pn,2n
[p_s_2n] = spectralExpansion2n(f_v_r,1);
%% Calculate Power Spectrum Pn,2n-1
f_rv_r_negative = rarImageIntensityCovariance(k,k_y,E_k,metadata,th,-1.*r);
[p_s_2n_1] = spectralExpansion2n_1(f_rv_r,f_rv_r_negative,f_v_r,1);
%% Calculate Power Spectrum Pn,2n-2
f_rv_r_zero = rarImageIntensityCovariance(k,k_y,E_k,metadata,th,0.*r);
[p_s_2n_2] = spectralExpansion2n_2(f_r_r,f_v_r,f_rv_r,f_rv_r_zero,f_rv_r_negative,1);
%% Full spectral expansion
%P_s = p_s_coeff.*(k_x.*beta).^m.*p;
P_s = zeros(2001,2001);
func = helperFunctions;
k_x_resize = func.resize(k_x,zeros(1,2001));
k_y_resize = func.resize(k_y,zeros(1,2001));
k_resize = func.resize(k,zeros(1,2001));
theta_resize = func.resize(theta,zeros(1,2001));
% for n=1:6
n=1;
    for m = 2*n-2:2*n
        fprintf('m = %d, ',m)
        fprintf('n = %d \n',n)
        coeff = ((k_x_resize.*beta).^m);
        if (m == 2*n)
            fprintf('P_(n,2n) used \n')
            P_s = P_s + coeff.*p_s_2n;
        end
        if (m == 2*n-1)
            fprintf('P_(n,2n-1) used \n')
            P_s = P_s + coeff.*p_s_2n_1;
        end        
        if (m == 2*n-2)
            fprintf('P_(n,2n-2) used\n')
            P_s = P_s + coeff.*p_s_2n_2;
        end
        
    end
% end
P_s = p_s_coeff.*P_s;
%% Plotting P_s
% Define for plotting
[t,r] = meshgrid(theta_resize,k_resize);
[x,y] = pol2cart(t,r);
figure;
contourf(x,y,real(P_s));
%%
figure;
contourf(k_x_resize, k_y_resize, real(P_s), 'LineStyle', 'none');
xlabel('k_x');
ylabel('k_y');
title('Spectrum P(k)');
colorbar;
%% Spectral expansions plots
func = helperFunctions;
k_x_resize = func.resize(k_x,zeros(1,2001));
k_y_resize = func.resize(k_y,zeros(1,2001));

figure;
plot(k_x_resize,abs(p_s_2n_1));
hold on;
plot(k_x_resize,abs(p_s_2n));
plot(k_x_resize,abs(p_s_2n_2));
hold off;
%colorbar;
% hold on;
% for i=1:5
%     [p_s_2n] = spectralExpansion2n(f_v_r,i);
%     display_name = [num2str(i), ': non linearity order'];
%     plot(k_x_resize,abs(p_s_2n),'DisplayName',display_name);
% end
% hold off;



%% Get sea surface elevation
% CDS data - 4 Aug 2022
seaLevelDataPath = "C:\Users\ryanj\Downloads\dataset-satellite-sea-level-global-49017144-e287-440f-a59d-61e7b6bf5c91\dt_global_twosat_phy_l4_20220804_vDT2021.nc";
sealevel = ncinfo(seaLevelDataPath);
sla = ncread(seaLevelDataPath,'sla');
lat_sla = ncread(seaLevelDataPath,'latitude');
long_sla = ncread(seaLevelDataPath,'longitude');
time_sla = ncread(seaLevelDataPath,'time');
time_sla = datetime(1950,1,1) + days(time_sla)
time_sla.Format = 'yyyy-MM-dd';
%% NOAA (Coastwatch) - 27 Sept 2023 data (Best so far) - Need DTU15

noaa_seaLevelPath = downloadNOAAWaveFile("https://coastwatch.noaa.gov/pub/socd/lsa/rads/sla/daily/nrt/2023/rads_global_nrt_sla_20230927_20230928_001.nc",'download.nc');
noaa_seaLevel = ncinfo(noaa_seaLevelPath);
sla = ncread(noaa_seaLevelPath,'sla');
lat_sla = ncread(noaa_seaLevelPath,'latitude');
long_sla = ncread(noaa_seaLevelPath,'longitude');
time_sla = ncread(noaa_seaLevelPath,'time');
time_sla = datetime(1950,1,1) + days(time_sla)
time_sla.Format = 'yyyy-MM-dd';
%% NOAA (S3A)
noaa_seaLevelPath = downloadNOAAWaveFile("https://www.star.nesdis.noaa.gov/data/pub0010/lsa/johnk/coastwatch/sa/sa_20230928.nc", 'download.nc');
noaa_seaLevel = ncinfo(noaa_seaLevelPath);
sla = ncread(noaa_seaLevelPath,'sla');
lat_sla = ncread(noaa_seaLevelPath,'lat');
long_sla = ncread(noaa_seaLevelPath,'lon');
time_sla = ncread(noaa_seaLevelPath,'time');
time_sla = datetime(1985,1,1) + seconds(time_sla)
time_sla.Format = 'yyyy-MM-dd';

% sla = sla(830:1876);
% lat_sla = lat_sla(830:1876);
% long_sla = long_sla(1:1876);
% time_sla = time_sla(830:1876);

%% Jason3 - https://podaac.jpl.nasa.gov/dataset/MERGED_TP_J1_OSTM_OST_CYCLES_V51?ids=&values=&search=measures%20v%205.1&provider=POCLOUD
noaa_seaLevelPath = downloadNOAAWaveFile("https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/MERGED_TP_J1_OSTM_OST_CYCLES_V51/Merged_TOPEX_Jason_OSTM_Jason-3_Cycle_1087.V5_1.nc", 'download.nc');
noaa_seaLevel = ncinfo(noaa_seaLevelPath);
ssha = ncread(noaa_seaLevelPath,'ssha'); % sourced in mm
mssh = ncread(noaa_seaLevelPath,'mssh'); % sourced in mm
sla = mssh + ssha;
lat_sla = ncread(noaa_seaLevelPath,'lat');
long_sla = ncread(noaa_seaLevelPath,'lon');
time_sla = ncread(noaa_seaLevelPath,'time');
time_sla = datetime(1992,1,1) + seconds(time_sla);
time_sla.Format = 'yyyy-MM-dd';

%% Get unique lat and long values
uniqueLat = unique(round(lat_sla,1));
uniqueLong = unique(round(long_sla,1));

% Initialize the sla matrix with NaN values
slaMatrix = nan(length(uniqueLat), length(uniqueLong));

% Assign sla values to the corresponding grid points
for i = 1:length(lat_sla)
    latIndex = find(uniqueLat == round(lat_sla(i),1));
    longIndex = find(uniqueLong == round(long_sla(i)),1);
    slaMatrix(latIndex, longIndex) = sla(i);
end

% Now you have slaMatrix as a matrix with corresponding sla values
%% Find sla val at desired lat and long
% Define the latitude and longitude coordinates where you want to find SLA
desiredLat = -34; % Replace with your desired latitude
desiredLong = 17; % Replace with your desired longitude

% Use interp2 to interpolate the SLA value
slaValue = interp2(uniqueLong, uniqueLat, slaMatrix, desiredLong, desiredLat);

% Check if the interpolated value is NaN (indicating no data at that point)
if isnan(slaValue)
    disp('No data available at the specified coordinates.');
else
    fprintf('SLA value at Lat %.3f, Long %.3f: %.3f mm\n', desiredLat, desiredLong, slaValue);
end

%%
figure(1)

m_proj('miller','lat',[min(uniqueLat(:)) max(uniqueLat(:))],...
'lon',[min(uniqueLong(:)) max(uniqueLong(:))])

% Next, plot the field using the M_MAP version of pcolor.

m_pcolor(uniqueLong,uniqueLat,slaMatrix);
shading flat;

% Add a coastline and axis values.

m_coast('patch',[.7 .7 .7])
m_grid('box','fancy')



% Add a colorbar and title.

colorbar
title('Sea Level Anomaly Height from NOAA (CoastWatch)');

%% Testing integrals
%dk = k(2) - k(1);
dth = theta(2) - theta(1);
M = [1 2];
int_th = cumsum(D)*dth;
int_th_trapz = cumtrapz(D)*dth;
figure;
plot(theta,int_th,'DisplayName', 'cumsum')
hold on;
plot(theta,int_th_trapz,'DisplayName' ,'cumtrapz');

