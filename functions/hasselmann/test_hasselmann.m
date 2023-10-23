%% Required files
% Import data
filepath = "D:\UCT\EEE4022S\Data\CPT\smaller_subset_incidence.nc";
%[Filename,filepath]=uigetfile('*.nc','Choose the exported .nc file from SNAP');
%vfilen=[filepath,Filename];
% Import data values
ncImport = ncinfo(filepath);
%VV_nc = ncread(filepath,'Sigma0_VV');
VV_nc = transectData_nc(:,:,3);
%% Required Metadata
metadata = ncinfo(filepath,'metadata');
th = ncread(filepath,'Incidence_Angle');
r = ones(512);
%% VV_nc redefine using transects and get correct th val
VV_nc = transectData_nc(:,:,1);
th = th(startPos_nc(1,1):startPos_nc(1,3),startPos_nc(1,2):startPos_nc(1,4));
%% Testing orbital velocity covariance function
th = incidenceAngle;
[f_v_r] = orbitalVelocityCovariance(k,k_y,E_k,metadata,r,th);
figure;
contour(k_x,k_y, 20*log10(abs(f_v_r)));
%surf(k_x,k_y, 20*log10(abs(f_v_r)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(f^v(r)))')
title('Orbital Velocity Covariance Function')
%% Testing autocovariance function of RAR image intensity
[f_r_r] = rarImageIntensityAutocovariance(k,k_y,E_k,metadata,r,th);
figure;
%contourf(k_x,k_y, 20*log10(abs(f_r_r)));
surf(k_x,k_y, 20*log10(abs(f_r_r)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(f^R(r)))')
title('RAR Image Intensity Autocovariance Function')
%% Testing covariance function of RAR image intensity and NL velocity
[f_rv_r] = rarImageIntensityCovariance(k,k_y,E_k,metadata,r,th);
figure;
%contourf(k_x,k_y, 20*log10(abs(f_rv_r)));
surf(k_x,k_y, 20*log10(abs(f_rv_r)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(f^{Rv}(r)))')
title('RAR image intensity and NL velocity Covariance Function')
%% Coefficient calculation
max_k = max(k);
min_k = min(k);
[p_s_coeff,beta,xi_sqr,k_x] = quasilinearCoeff(k,k_y,k_x,E_k,metadata,th);
%%
figure;
%contour(k_x,k_y, 20*log10(abs(p_s_coeff)));
%surf(k_x,k_y, 20*log10(abs(p_s_coeff)),'lineStyle','none');
plot(k_x,p_s_coeff);
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(exp(-k_x^2\xi^2))')
title('exp(-k_x^2\xi^2)')
%% Calculate Power Spectrum Pn,2n
[p_s_2n] = spectralExpansion2n(f_v_r,1);
figure;
%contour(k_x,k_y, 20*log10(abs(p_s_2n)));
surf(k_x,k_y, 20*log10(abs(p_s_2n)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(P^S_{n,2n}))')
title('Power Spectrum, P^S_{n,2n}')
%% Calculate Power Spectrum Pn,2n-1
f_rv_r_negative = rarImageIntensityCovariance(k,k_y,E_k,metadata,-1.*r,th);
[p_s_2n_1] = spectralExpansion2n_1(f_rv_r,f_rv_r_negative,f_v_r,1);
figure;
contour(k_x,k_y, 20*log10(abs(p_s_2n_1)));
%surf(k_x,k_y, 20*log10(abs(p_s_2n_1)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
%zlabel('20log(abs(P^S_{n,2n-1}))')
title('Power Spectrum, P^S_{n,2n-1}')
%% Calculate Power Spectrum Pn,2n-2
f_rv_r_zero = rarImageIntensityCovariance(k,k_y,E_k,metadata,0.*r,th);
[p_s_2n_2] = spectralExpansion2n_2(f_r_r,f_v_r,f_rv_r,f_rv_r_zero,f_rv_r_negative,1);
figure;
contour(k_x,k_y, 20*log10(abs(p_s_2n_2)));
%surf(k_x,k_y, 20*log10(abs(p_s_2n_2)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
%zlabel('20log(abs(P^S_{n,2n-2}))')
title('Power Spectrum, P^S_{n,2n-2}')
%colorbar;
%% Full spectral expansion
%P_s = p_s_coeff.*(k_x.*beta).^m.*p;
P_s = zeros(512,512);
func = helperFunctions;
k_x_resize = func.resize(k_x,zeros(1,2001));
k_y_resize = func.resize(k_y,zeros(1,2001));
k_resize = sqrt(k_x_resize.^2+k_y_resize.^2);
theta_resize = func.resize(theta',zeros(1,2001));
%for n=1:6
n=1;
    for m = 2*n-2:2*n
        fprintf('m = %d, ',m)
        fprintf('n = %d \n',n)
        coeff = ((k_x.*beta).^m);
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
%end
P_s = p_s_coeff.*P_s;

figure;
%contour(k_x,k_y, 20*log10(abs(P_s)));
surf(k_x,k_y, 20*log10(abs(P_s)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
%zlabel('20log(abs(P^S))')
title('Power Spectrum, P^S')
%colorbar;
%% Testing filter functions - HRK
xi_sqr = 70.^2;
mu_Th = 0.5;
look = func.getLook(metadata);
look = func.look(look);
polarisation = func.getPolarisation(metadata);
th = incidenceAngle;

% Need to check resizing is correct
%k_new = func.resize(k,th(1,:));
k_new = k;
%k_y = func.resize(k_y,th(1,:));
k_l = func.kl(look,k_y);
k_l_inv = func.kl(look,k_y_inv);
omega = func.omega(k_new);

Th_k = func.hydroMTF(omega,mu_Th,k_new,k_y);
Tt_k = func.tiltMTF(polarisation,k_l,th);
TR_k = func.rarMTF(Tt_k,Th_k);

Th_k_inv = func.hydroMTF(omega,mu_Th,k_inv,k_y_inv);
Tt_k_inv = func.tiltMTF(polarisation,k_l_inv,th);
TR_k_inv = func.rarMTF(Tt_k_inv,Th_k_inv);

HRK = exp(-k_x.*xi_sqr).*TR_k;
HRK = func.resize(HRK(:,2:end),th);
HRK = fftshift(fft(HRK));
figure;
%surf(10*log10(abs(HRK)),'LineStyle','none');
contour(20*log10(abs(HRK)));
%imagesc(abs(HRK))

%% HVBk

Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tvb_k = func.velocityBunchingMTF(beta,k_x,Tv_k);
Tv_k_inv = func.rangeVelocityTF(omega,th,k_l,-k_new);
Tvb_k_inv = func.velocityBunchingMTF(beta,k_x,Tv_k_inv);

HVbK = exp(-k_x.*xi_sqr).*abs(Tvb_k).^2;
HVbK = func.resize(HVbK(:,2:end),th);
HVbK = fftshift(fft(HVbK));
figure;
%surf(10*log10(abs(HVbK)),'LineStyle','none');
contour(10*log10(abs(HVbK)));
%imagesc(HVbK)

%% HSK
Ts_k_inv = func.sarImagingMTF(TR_k_inv,Tvb_k_inv);
Ts_k = func.sarImagingMTF(TR_k,Tvb_k);
HSK = exp(-k_x.*xi_sqr).*abs(Ts_k).^2; 
HSK = fftshift(fft(HSK));
figure;
%surf(10*log10(abs(HSK)),'LineStyle','none');
contour(10*log10(abs(HSK)));
%imagesc(HSK)

%% HintK
HintK = exp(-k_x.*xi_sqr).*(TR_k.*conj(Tvb_k) + conj(TR_k).*Tvb_k); 
HintK = fftshift(fft(HintK));
figure;
%surf(10*log10(abs(HintK)),'LineStyle','none');
contour(10*log10(abs(HintK)));
%imagesc(HintK)
%% Inversion
mu = calculateWeight(intensityFFT);
B = calculateNormalisationConstant(E_k);
%J = costFunction(P_s_pipeline,sarData,E_k,B,mu);
%% Plotting P_s
% Define for plotting
[t,r] = meshgrid(theta,k);
[x,y] = pol2cart(t,r);
figure;
%plot(t,r)
contour(x,y,20*log10(abs(P_s)));
yline(0);
xline(0);
grid on;
xlabel('\theta (rad)'), ylabel('\theta (rad)');
%%
figure;
%surf(k_x, k_y, 20*log10(abs(P_s)),'LineStyle', 'none');
contour(k_x, k_y, 20*log10(abs(P_s)));
xlabel('k_x');
ylabel('k_y');
grid on;
title('Spectrum Abs(P(k))');
%xticks(-10:0.1:10)
%set(gca,'XTick',-10:0.1:10) 
%set(gca,'YTick',-10:0.1:10)
%colorbar;
%%
mean_VV = mean(VV_nc);
cvar = var((VV_nc-mean_VV)./mean_VV);
%%
figure;
imagesc(20*log10(abs(P_s)));
%imagesc(VV_nc);
colorbar;
%% Check linear SAR spectrum
func = helperFunctions;

% Tv_k = func.rangeVelocityTF(omega,incidenceAngle,k_l,k);
% Tvb_k = func.velocityBunchingMTF(beta,k_x,Tv_k);
% TS_k = func.sarImagingMTF(TR_k,Tvb_k);

P_s_lin = imageVarianceSpectrum(k,k_x,k_y,k_inv,k_x_inv,k_y_inv,E_k,E_k_inv,metadata,th);
figure;
%contour(k_x,k_y, 20*log(abs(P_s_lin)));
surf(k_x,k_y, 20*log10(abs(P_s_lin)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(P^S_k))')
title('Image Variance Spectrum')

%% Get SAR intensity power spectrum
% Calculate the 2D Fourier transform of the intensity data
intensityFFT = abs(fftshift(fft2(VV_nc)));
%%
% Plot the power spectrum
dk_x = 0.005890997229533;
dk_y = -0.005226593472072;
above0intesity = intensityFFT;
above0intesity(intensityFFT < 1) = 0;
%%
figure;
%contourf(k_x./dk_x,k_y./dk_y,20*log10(intensityFFT))
contourf(20*log10(intensityFFT))
%surf(20*log10(intensityFFT),'LineStyle','none')
colorbar;
%imagesc(20*log10(intensityFFT))
%imshow(VV_nc)
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

