%clear
clc;
%% Determine wave file to download
func = helperFunctions;
captureDate = func.getCaptureDate(metadata);
[noaaDateStr, noaaHourStr] = func.getNOAAParams(captureDate);
noaaUrl = ['https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.',noaaDateStr,'/',noaaHourStr,'/wave/gridded/gdaswave.t',noaaHourStr,'z.gsouth.0p25.f000.grib2'];
%%
% Hardcode as SAR data not available in window
noaaUrl = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.20231010/12/wave/gridded/gdaswave.t12z.gsouth.0p25.f000.grib2';
filePath = downloadNOAAWaveFile(noaaUrl,'wave_data.grib2');
%% Real Deal
% NOAA SOUTH IS NOT DOWNLOADING SOUTH - USE DOWNLOAD IN repo NOT functions
% Can only have 1 .grib2 file in path
outStruct = getGribStruct("C:\Users\ryanj\Downloads","C:\Users\ryanj\OneDrive - University of Cape Town\4. Fourth Year\Second Semester\EEE4022S\repo\sar-parameter-extraction-pipeline\functions\");

%% Get variables
% Define constants
g = 9.81;

latStart = grid_lat(1,3)
latEnd = grid_lat(1,3)
lonStart = grid_lon(1,1)
lonEnd = grid_lon(1,2)
%%
% Get range of values in lat and long
waveVals = getSubsetWaveVals(outStruct,latStart,lonStart,latEnd,lonEnd);

%% Calculate Wave Spectrum Params
% lambda = g.*T.^2/2*pi;
% k = (4*pi^2)./(g.*testWavePeriod(1,1)); % Check if correct
% k_x = k.*cos(testWaveDirection(1,1));
% k_y = k.*sin(testWaveDirection(1,1));

%% Plot Wave Data on world map
noaaDataPlot('miller',outStruct,'windDirection')
%% 3D plot
% Instantiate variables
image_size = 512;
%image_size = size(VV_nc,1);
%w = (0:0.025:2*pi)';
w = linspace(0,2*pi,image_size)';
f = linspace(0,1,image_size)';
Hs = waveVals.significantWaveHeight(1,1);
T0 = waveVals.significantWavePeriod(1,1);
w0 = 2*pi./T0;
f0 = 1./T0;
%% Calculate S(\omega)
%gamma_val = 1.308.*ones(size(Hs));
gamma_val = 1.308;
%%
S = generateSingleJONSWAP(Hs,w0,gamma_val,w);
% gamma_val = 1.308;
% sigma = 0.07; % or 0.09
% nonLin = exp(3.3-0.03.*(gamma_val-1)./sigma)
%S = wavespec(7,[Hs,w0,gamma_val],w,0);
figure;
plot(w,S);
%plot(f,S);
%title(['One-dimensional wave spectrum, E(\omega) at ', num2str(waveVals.latitude(1,1)), 'S, ', num2str(waveVals.longitude(1,1)), 'E'])
legend('JONSWAP for \gamma = 1.308','Location','Northeast')
xlabel('\omega [rad/s]')
ylabel('E(\omega) [m^2/rad/Hz]')
grid on;
% set(gca,'XTick',0:0.25:1) 
% set(gca,'XTickLabel',{'0','0.25','0.5','0.75','1'})
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
set(gca,'FontSize',12)
legend('show')
%matlab2tikz('../plots/oneDWaveSpec.tex');
%% Calcuate S(\omega) for multiple lats and longs
S = generateMultipleJONSWAP(waveVals,gamma_val,w,1);
%% Sspectral bandwidth calculation
S_1 = S(:,:,1);
S_1 = S_1(2:end);
S_1 = func.resize(S_1(2:end),w);
S_norm = S_1 / (trapz(S_1).*dw);
spectral_mean = trapz(w.*S_norm).*dw;
spectral_var = trapz((w-spectral_mean).^2.*S_norm).*dw;
spectral_bw = sqrt(spectral_var);

%% 2D plot of S(\omega)
figure;
%plot(w,S/(T0*Hs^2))
%title('S(\omega)/(H_s^2 T_0)')

hold on;
m = length(waveVals.longitude);
for i = 1:2
    lat_index = ceil(i/m);
    lon_index = mod(i-1,m)+1;
    display_name = [num2str(waveVals.latitude(lat_index)), 'S, ' num2str(waveVals.longitude(lon_index)), 'E'];
    % For plotting multiple lat vals
    % switch lat_index
    %     case 1
    %         colourToPlot = "#0072BD";
    %     case 2
    %         colourToPlot = "#EDB120";
    %     case 3 
    %         colourToPlot = "#A2142F";
    % end
    % plot(w,S(:,i),'DisplayName',display_name, 'Color',colourToPlot);
    plot(w,S(:,i),'DisplayName',display_name);
    maxIndex = find(S(:,i) == max(S(:,i)))
    switch lon_index
        case 1
            colourToPlot = "#0072BD";
        case 2              
            colourToPlot = "#D95319";
    end
    S(maxIndex,i)
    xline(w(maxIndex),LineWidth=1,DisplayName=['\omega = ',num2str(round(w(maxIndex),3))],Color=colourToPlot,LineStyle="--");
    %legend(["Wave spectrum at: ", num2str(testLat(lat_index)), ", " num2str(testLon(lon_index))]);
    %fprintf('Iteration %d: Latitude %.2f, Longitude %.2f\n', i, testLat(lat_index), testLon(lon_index));
end
hold off

%figure;
%plot(w,S(:,1));
%plot(f,S);
%title('E(\omega)')
%legend('show');
%legend('JONSWAP for \gamma = 1.308','Location','Northeast')
xlabel('\omega [rad/s]')
ylabel('E(\omega) [m^2/rad/Hz]')
grid on;
% set(gca,'XTick',0:0.25:1) 
% set(gca,'XTickLabel',{'0','0.25','0.5','0.75','1'})
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
set(gca,'FontSize',12)
legend('show')
%matlab2tikz('../plots/valWaveSpec_toShore_1.tex');
%% NonLinearity play around
f_m = f0
nLinOrder = (log(10).*(2*pi.*f_m)^4)./gamma_val.^2


%% Directional spreading function, D(\omega)
%update to take in desired lat and long coords
[D,theta] = generateDirectionalDistribution(waveVals,w,1);
E = generate2DWaveSpectrum(S,D);
%%
Tsig = waveVals.significantWavePeriod(1,1);
Tp = Tsig./0.95;
wsig = 2*pi./Tsig;
wp = 2*pi./Tp;
sig_th = 26.9.*((wsig)/(wp)).^(-1.05); % in degrees
for i=1:length(w)
    sig_th(i) = 26.9.*((wsig)/(wp)).^(-1.05); % in degrees
    if (w(i)>=wp)
        sig_th(i) = 26.9.*((wsig)/(wp)).^(0.68);
    end
end
%sig_th = 30;
s = 2./(deg2rad(sig_th).^2) - 1;
%s = 84;
Tp = waveVals.significantWavePeriod(1,1)./0.95;
%m = 0.0585.*Tp.^2 - 0.3988.*Tp + 3.0546;
%theta = (-pi:0.025:pi)'; % Old definition of theta
th_wave = waveVals.direction(1,1);
th_wind = waveVals.windDirection(1,1);

low = th_wave - pi;
high = th_wave + pi;
L = low  + rem(th_wave - low,0.025);
H = high - rem(high - th_wave,0.025);
theta = linspace(L,H,252)';

%A1 = gamma(0.5*m+1)./(gamma(0.5*m+1/2).*sqrt(pi));
%D = A1.*cos(theta).^m;
i = 1;
j = 1;
while (i<=252)
    A2 = gamma(s(i)+1)./(gamma(s(i)+1/2).*2.*sqrt(pi));
    D(:,j) = abs(A2 .* cos(0.5.*(th_wave-theta)).^(2*s(i))); % abs because D is complex. 
    %t = ones(size(S)); 
    E(:,:,:,j) = S .* D(:,j)';
    i = i + 251;
    j = j + 1;
end

%A2 = gamma(s+1)./(gamma(s+1/2).*2.*sqrt(pi));
%D = abs(A2 .* cos(0.5.*(th_wave-theta)).^(2*s)); % abs because D is complex. 
%t = ones(size(S)); 
%E = S .* D';
%%
% Replace values close to zero with NaN (e.g., values smaller than a threshold)
threshold = 1e-4;  % Adjust the threshold as needed

% Create a copy of the original matrix and apply the threshold
E_cleaned = E;
E_cleaned(abs(E_cleaned) < threshold) = NaN;

% Define for plotting
[t,r] = meshgrid(theta,w);
[x,y] = pol2cart(t,r);



figure;
i = 1;
j = 1;
hold on;
for j=1:2
    display_name = ['s = ', num2str(s(i))];
    plot(theta,D(:,j),'DisplayName',display_name);
    i = i + 511;
    j+1;
end
%plot(theta,D);
%title('Directional Distribution Function');
%title(['D(\theta) for s = ', num2str(s)]);
xlabel('\theta [rad]');
ylabel('D(\theta) [1/rad]')
grid on;
set(gca,'defaultAxesTickLabelInterpreter','latex'); 
set(gca,'XTick',L:pi/2:H) 
set(gca,'XTickLabel',{'','\theta_{wave}-\pi/2','\theta_{wave}','\theta_{wave}+\pi/2',''})
set(gca,'FontSize',12)
legend('show')
%matlab2tikz('../plots/directionalSpreading.tex');
%% Contour of 2D wave spectrum
E = E(:,:,1,2);
%%
figure;
contour(x,y,E)
yline(0);
xline(0);
grid on;
xlabel('\omega [rad/s]'), ylabel('\omega [rad/s]');
set(gca,'XTick',-2*pi:pi/2:2*pi) 
set(gca,'XTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
set(gca,'YTick',-2*pi:pi/2:2*pi) 
set(gca,'YTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
set(gca,'FontSize',12)
%matlab2tikz('../plots/contour_E_omega.tex');
%% 3D Spectrum
figure;
h=surf(x,y,E);
grid on;
xlabel('\omega [rad/s]'), ylabel('\omega [rad/s]');
set(gca,'XTick',-2*pi:pi/2:2*pi) 
set(gca,'XTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
set(gca,'YTick',-2*pi:pi/2:2*pi) 
set(gca,'YTickLabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi', '3\pi/2', '2\pi'})
set(gca,'FontSize',22)
set(h,'LineStyle','none');
c = colorbar;
%c.Label.String = barStr;
hL = ylabel(c,'[m^2/rad/Hz]');
set(hL,'Rotation',0);
% figure(6);
% h = surf(w,theta,E_cleaned);
% set(h,'LineStyle','none');
% xlabel('\omega (rad/s)')
% ylabel('\theta (rad)')
zlabel('E(\omega,\theta)')
%matlab2tikz('../plots/surf_E_omega.tex');

%% Calculate E(k_x,k_y)
d = 70;
k = (w.^2./g)';
[E_k,k] = waveNumberSpectrum(E,w,k,d);
[E_k_inv,k_inv] = waveNumberSpectrum(E,w,-k,d);
%%
d = 70; % 70m depth at CP
%w = sqrt(g.*k.*tanh(k.*d)); % (eq. 5.4.17 in holthuisjen)
c = sqrt((g./k).*tanh(k.*d)); % tanh in radians (eq. 5.4.23 in holthuisjen) output in m/s^2
n = 0.5*(1+(2.*k.*d)./sinh(2.*k.*d)); % sinh in radians (eq. 5.4.32 in holthuisjen)
c_g = n.*c;
test = (c.*c_g)./w;
E_k = ((c.*c_g)./w).*E;

%% k definition of E(k_x,k_y)
%k = linspace(0,10,252);
%k = (w.^2./g)';
k_x = k.*cos(waveVals.direction(2,2));
k_y = k.*sin(waveVals.direction(2,2));

k_x_inv = k_inv.*cos(waveVals.direction(2,2));
k_y_inv = k_inv.*sin(waveVals.direction(2,2));

%% Validation of n
nValidation(k,70,5);

%% Plots of E(k_x,k_y)

% Contour of 2D wave spectrum
%E = E(:,:,2);
figure;
contour(k_x,k_y,E_k)
grid on;
yline(0);
xline(0);
xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex');
%matlab2tikz('../plots/co
%set(gca,'XTick',-pi/2:pi/6:pi/2) 
%set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
%set(gca,'YTick',-pi/2:pi/6:pi/2) 
%set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
set(gca,'FontSize',12)
%matlab2tikz('../plots/Valcontour_E_k_inv.tex');
%% 3D Spectrum
figure;
h=surf(k_x,k_y,E_k);
grid on;
xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex'),zlabel('$E(k_{x},k_{y})$','interpreter','latex');
set(h,'LineStyle','none');
c = colorbar;
%c.Label.String = barStr;
hL = ylabel(c,'[m^2/rad/Hz]');
set(hL,'Rotation',0);
set(gca,'FontSize',32)
%set(gca,'XTick',-pi/2:pi/6:pi/2) 
%set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
%set(gca,'YTick',-pi/2:pi/6:pi/2) 
%set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
%matlab2tikz('../plots/surf_E_k.tex');
