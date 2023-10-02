%clear
clc;

filePath = downloadNOAAWaveFile('https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gdas.20230923/12/wave/gridded/gdaswave.t12z.gsouth.0p25.f000.grib2','wave_data.grib2');
%% Real Deal
% NOAA SOUTH IS NOT DOWNLOADING SOUTH - USE DOWNLOAD IN repo NOT functions
% Can only have 1 .grib2 file in path
outStruct = getGribStruct("C:\Users\ryanj\Downloads","C:\Users\ryanj\OneDrive - University of Cape Town\4. Fourth Year\Second Semester\EEE4022S\repo\sar-parameter-extraction-pipeline\functions\");

%% Get variables
% Define constants
g = 9.81;

% Get range of values in lat and long
waveVals = getSubsetWaveVals(outStruct,-34,17.25,-34,17.75);

%% Calculate Wave Spectrum Params
% lambda = g.*T.^2/2*pi;
% k = (4*pi^2)./(g.*testWavePeriod(1,1)); % Check if correct
% k_x = k.*cos(testWaveDirection(1,1));
% k_y = k.*sin(testWaveDirection(1,1));

%% Plot

noaaDataPlot('miller',outStruct,'windSpeed')
%%
figure(1)

m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
'lon',[min(lon(:)) max(lon(:))])

% Next, plot the field using the M_MAP version of pcolor.

m_pcolor(lon,lat,waveHeight.');
shading flat;

% Add a coastline and axis values.

m_coast('patch',[.7 .7 .7])
m_grid('box','fancy')



% Add a colorbar and title.

colorbar
title('WAVEWATCH III Wave Height from NOAA (NCEP)');
%% 3D plot
% Instantiate variables
w = (0:0.025:2*pi)';
f = linspace(0,1,252)';
Hs = testWaveHeight(1,1);
T0 = testWavePeriod(1,1);
w0 = 2*pi./T0;
f0 = 1./T0;
%% Calculate S(\omega)
%gamma_val = 1.308.*ones(size(Hs));
gamma_val = 1.308;
sigma = 0.07; % or 0.09
nonLin = exp(3.3-0.03.*(gamma_val-1)./sigma)
S = wavespec(7,[Hs,w0,gamma_val],w,0);
%% Calcuate S(\omega) for multiple lats and longs
w = (0:0.025:2*pi)';
Hs = testWaveHeight;
T0 = testWavePeriod;
w0 = 2*pi./T0;
f0 = 1./T0;
gamma_val = 1.308;
count = 1;
for i = 1:length(testLat)
    for j = 1:length(testLon)
        Hs = testWaveHeight(j,i);
        T0 = testWavePeriod(j,i);
        w0 = 2*pi./T0;
        f0 = 1./T0;
        S(:,:,count) = wavespec(7,[Hs,w0,gamma_val],w,0);
        count = count + 1;
    end
end
%% 2D plot of S(\omega)
figure(2);
%plot(w,S/(T0*Hs^2))
%title('S(\omega)/(H_s^2 T_0)')
hold on;
m = length(testLon);
for i = 1:count-1
    lat_index = ceil(i/m);
    lon_index = mod(i-1,m)+1;
    display_name = [num2str(testLat(lat_index)), 'S, ' num2str(testLon(lon_index)), 'E'];
    plot(w,S(:,i),'DisplayName',display_name);
    %legend(["Wave spectrum at: ", num2str(testLat(lat_index)), ", " num2str(testLon(lon_index))]);
    %fprintf('Iteration %d: Latitude %.2f, Longitude %.2f\n', i, testLat(lat_index), testLon(lon_index));
end
hold off

%figure;
%plot(w,S(:,1));
%plot(f,S);
title('S(\omega)')
%legend('show');
%legend('JONSWAP for \gamma = 1.308','Location','Northeast')
xlabel('\omega (rad/s)')
grid on;
% set(gca,'XTick',0:0.25:1) 
% set(gca,'XTickLabel',{'0','0.25','0.5','0.75','1'})
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
set(gca,'FontSize',12)
legend('show')
%matlab2tikz('../plots/waveSpec_toShore.tex');
%% NonLinearity play around
f_m = f0
nLinOrder = (log(10).*(2*pi.*f_m)^4)./gamma_val.^2


%% Directional spreading function, D(\omega)
Tsig = testWavePeriod(1,1);
Tp = Tsig./0.95;
sig_th = 26.9.*((1./Tsig)/(1./Tp)).^(-1.05); % in degrees
%sig_th = 30;
s = 2./(deg2rad(sig_th).^2) - 1;
%s = 84;
Tp = testWavePeriod(1,1)./0.95;
m = 0.0585.*Tp.^2 - 0.3988.*Tp + 3.0546;
%theta = (-pi:0.025:pi)'; % Old definition of theta
th_wave = testWaveDirection(1,1);
th_wind = testWindDirection(1,1);

low = th_wave - pi;
high = th_wave + pi;
L = low  + rem(th_wave - low,0.025);
H = high - rem(high - th_wave,0.025);
theta = linspace(L,H,252)';

%A1 = gamma(0.5*m+1)./(gamma(0.5*m+1/2).*sqrt(pi));
%D = A1.*cos(theta).^m;
A2 = gamma(s+1)./(gamma(s+1/2).*2.*sqrt(pi));
D = abs(A2 .* cos(0.5.*(th_wave-theta)).^(2*s)); % abs because D is complex. 
%t = ones(size(S)); 
E = S .* D';
% Replace values close to zero with NaN (e.g., values smaller than a threshold)
threshold = 1e-4;  % Adjust the threshold as needed

% Create a copy of the original matrix and apply the threshold
E_cleaned = E;
E_cleaned(abs(E_cleaned) < threshold) = NaN;

% Define for plotting
[t,r] = meshgrid(theta,w);
[x,y] = pol2cart(t,r);



figure(3);
plot(theta,D);
title(['D(\theta) for s = ', num2str(s)]);
xlabel('\theta (rad)');
ylabel('D(\theta) (1/rad)')
grid on;
set(gca,'defaultAxesTickLabelInterpreter','latex'); 
set(gca,'XTick',L:pi/2:H) 
set(gca,'XTickLabel',{'','\theta_{wave}-\pi/2','\theta_{wave}','\theta_{wave}+\pi/2',''})
set(gca,'FontSize',12)
%matlab2tikz('../plots/directionalSpreading.tex');
%% Contour of 2D wave spectrum
E = E(:,:,2);
figure(4);
contour(x,y,E)
grid on;
xlabel('\omega (rad/s)'), ylabel('\omega (rad/s)');
%set(gca,'XTick',-pi/2:pi/6:pi/2) 
%set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
%set(gca,'YTick',-pi/2:pi/6:pi/2) 
%set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
set(gca,'FontSize',12)
%matlab2tikz('../plots/contour_E_omega.tex');
% 3D Spectrum
figure;
h=surf(x,y,E);
grid on;
xlabel('\omega (rad/s)'), ylabel('\omega (rad/s)'),zlabel('E(\omega,\theta)');
set(h,'LineStyle','none');
%set(gca,'XTick',-pi/2:pi/6:pi/2) 
%set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
%set(gca,'YTick',-pi/2:pi/6:pi/2) 
%set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
colorbar;
%matlab2tikz('../plots/surf_E_omega.tex');
% figure(6);
% h = surf(w,theta,E_cleaned);
% set(h,'LineStyle','none');
% xlabel('\omega (rad/s)')
% ylabel('\theta (rad)')
% zlabel('E(\omega,\theta)')
% colorbar;

%% Calculate E(k_x,k_y)
d = 1; % 70m depth at CP
%w = sqrt(g.*k.*tanh(k.*d)); % (eq. 5.4.17 in holthuisjen)
c = sqrt((g./k).*tanh(k.*d)); % tanh in radians (eq. 5.4.23 in holthuisjen) output in m/s^2
n = 0.5*(1+(2.*k.*d)./sinh(2.*k.*d)); % sinh in radians (eq. 5.4.32 in holthuisjen)
c_g = n.*c;
test = (c.*c_g)./w;
E_k = ((c.*c_g)./w).*E;

%% k definition of E(k_x,k_y)
%k = linspace(0,10,252);
k = (w.^2./g)';
k_x = k.*cos(testWaveDirection(1,1));
k_y = k.*sin(testWaveDirection(1,1));

%% Validation of n
nValidation(k,70,5);

%% Plots of E(k_x,k_y)

% Contour of 2D wave spectrum
%E = E(:,:,2);
figure(6);
contour(k_x,k_y,E_k)
grid on;
xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex');
%set(gca,'XTick',-pi/2:pi/6:pi/2) 
%set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
%set(gca,'YTick',-pi/2:pi/6:pi/2) 
%set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
set(gca,'FontSize',12)
%matlab2tikz('../plots/contour_E_k.tex');
% 3D Spectrum
figure;
h=surf(k_x,k_y,E);
grid on;
xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex'),zlabel('$E(k_{x},k_{y})$','interpreter','latex');
set(h,'LineStyle','none');
%set(gca,'XTick',-pi/2:pi/6:pi/2) 
%set(gca,'XTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
%set(gca,'YTick',-pi/2:pi/6:pi/2) 
%set(gca,'YTickLabel',{'-\pi/2','-\pi/3','-\pi/6','0','\pi/6', '\pi/3', '\pi/2'})
colorbar;
%matlab2tikz('../plots/surf_E_k.tex');
