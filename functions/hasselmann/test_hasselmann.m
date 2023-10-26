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
%% 
r = ones(512);
%% VV_nc redefine using transects and get correct th val
VV_nc = transectData_nc(:,:,1);
th = th(startPos_nc(1,1):startPos_nc(1,3),startPos_nc(1,2):startPos_nc(1,4));
%% Testing orbital velocity covariance function
th = incidenceAngle(:,:,1);
[f_v_r] = orbitalVelocityCovariance(k,k_y,E_k,metadata,r,th);
SARSpectrumPlot(abs(f_v_r),0,spectralBW,1);
%%
figure;
%contour(k_x,k_y, 20*log10(abs(f_v_r)));
surf(k_x,k_y, 20*log10(abs(f_v_r)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(f^v(r)))')
title('Orbital Velocity Covariance Function')
%% Testing autocovariance function of RAR image intensity
[f_r_r] = rarImageIntensityAutocovariance(k,k_y,E_k,metadata,r,th);
SARSpectrumPlot(abs(f_r_r),0,spectralBW,1);
%%
figure;
%contourf(k_x,k_y, 20*log10(abs(f_r_r)));
surf(k_x,k_y, 20*log10(abs(f_r_r)),'lineStyle','none');
xlabel('k_x')
ylabel('k_y')
zlabel('20log(abs(f^R(r)))')
title('RAR Image Intensity Autocovariance Function')
%% Testing covariance function of RAR image intensity and NL velocity
[f_rv_r] = rarImageIntensityCovariance(k,k_y,E_k,metadata,r,th);
SARSpectrumPlot(abs(f_rv_r),0,spectralBW,1);
%%
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
[p_s_coeff,beta,xi_sqr,~] = quasilinearCoeff(k,k_y,k_x,E_k,metadata,th);
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
SARSpectrumPlot(abs(p_s_2n),120,spectralBW,1);
%contour(X,Y, 20*log10(abs(subSpectralExpansion2n)));
%surf(k_x,k_y, 20*log10(abs(p_s_2n)),'lineStyle','none');
% xlabel('k_x')
% ylabel('k_y')
% zlabel('20log(abs(P^S_{n,2n}))')
title('Power Spectrum, P^S_{n,2n}')
%matlab2tikz('../plots/Pn_2n.tex');
%% Calculate Power Spectrum Pn,2n-1
f_rv_r_negative = rarImageIntensityCovariance(k,k_y,E_k,metadata,-1.*r,th);
[p_s_2n_1] = spectralExpansion2n_1(f_rv_r,f_rv_r_negative,f_v_r,1);
SARSpectrumPlot(abs(p_s_2n_1),120,spectralBW,1);
title('Power Spectrum, P^S_{n,2n-1}')
%matlab2tikz('../plots/Pn_2n_1.tex');
%%
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
SARSpectrumPlot(abs(p_s_2n_2),0,spectralBW,1)
%matlab2tikz('../plots/Pn_2n_2_NoThreshold.tex');
%%
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
%theta_resize = func.resize(theta',zeros(1,2001));
%for n=1:6
n=1;
    for m = 2*n-2:2*n
        fprintf('m = %d, ',m)
        fprintf('n = %d \n',n)
        coeff = ((k_x.*beta).^m);
        if (m == 2*n)
            fprintf('P_(n,2n) used \n')
            % P_s = P_s + coeff.*p_s_2n;
            P_s = coeff.*p_s_2n;
            SARSpectrumPlot(abs(P_s),0,spectralBW,1);
            title(['m = ',num2str(m),', n = ',num2str(n)]);
            %matlab2tikz('../plots/Pn_2n.tex');
        end
        if (m == 2*n-1)
            fprintf('P_(n,2n-1) used \n')
            % P_s = P_s + coeff.*p_s_2n_1;
            P_s = coeff.*p_s_2n_1;
            SARSpectrumPlot(abs(P_s),0,spectralBW,1);
            % title(['m = ',num2str(m),', n = ',num2str(n)]);
            %matlab2tikz('../plots/Pn_2n_1.tex');
        end        
        if (m == 2*n-2)
            fprintf('P_(n,2n-2) used\n')
            P_s = P_s + coeff.*p_s_2n_2;
            % P_s = coeff.*p_s_2n_2;
            SARSpectrumPlot(abs(P_s),0,spectralBW,1);
            title(['m = ',num2str(m),', n = ',num2str(n)]);
            %matlab2tikz('../plots/Pn_2n_2_coeff.tex');
        end
        
    end
%end
P_s = p_s_coeff.*P_s;
%% 
coeff = ((k_x.*beta).^2);
k_x_pos = -k_x;
dkx = 0.002511482945409;
dky = 0.002511482945409;
reflected_k_x = [fliplr(k_x), -k_x];
reflected_k_y = [fliplr(k_y), -k_y]';
duplicate_k_y = [fliplr(k_y), k_y]';
reflected_coeff = [fliplr(coeff), coeff];
reflected_coeff = abs(reflected_coeff .*  duplicate_k_y);
figure;
%plot(reflected_k_x,reflected_coeff);
surf(reflected_k_x/dkx,reflected_k_y/dky,reflected_coeff,'LineStyle','none');
xlabel('$k_x/\Delta k$','In terpreter','latex');
ylabel('$k_y/\Delta k$','Interpreter','latex');
zlabel('$(k_{x} \cdot \beta)^2$','Interpreter','latex')
grid on;
grid minor;
set(gca,'FontSize',40)
%%
figure;
plot(k_x,coeff);
xlabel('$k_x$','Interpreter','latex');
ylabel('$(k_{x} \cdot \beta)^m$','Interpreter','latex')
grid on;
grid minor;
%matlab2tikz('../plots/powerSpectrum_coeff.tex');
%%
figure;
plot(k_x,coeff);
xlabel('$k_x$','Interpreter','latex');
ylabel('$(k_{x} \cdot \beta)^m$','Interpreter','latex')
grid on;
grid minor;
%matlab2tikz('../plots/powerSpectrum_coeff.tex');
%%
figure;
contour(k_x,k_y, 20*log10(abs(P_s)));
%surf(k_x,k_y, 20*log10(abs(P_s)),'lineStyle','none');
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
SARSpectrumPlot(abs(HRK),120,spectralBW,1);
title(['HRK Filter function']);

%% HVBk

Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tvb_k = func.velocityBunchingMTF(beta,k_x,Tv_k);
Tv_k_inv = func.rangeVelocityTF(omega,th,k_l,fliplr(-k_new));
Tvb_k_inv = func.velocityBunchingMTF(beta,k_x_inv,Tv_k_inv);

HVbK = exp(-k_x.*xi_sqr).*abs(Tvb_k).^2;
HVbK = func.resize(HVbK(:,2:end),th);
HVbK = fftshift(fft(HVbK));
SARSpectrumPlot(abs(HVbK),120,spectralBW,1);
title(['HVbK Filter function']);

%% HSK
Ts_k_inv = func.sarImagingMTF(TR_k_inv,Tvb_k_inv);
Ts_k = func.sarImagingMTF(TR_k,Tvb_k);
HSK = exp(-k_x.*xi_sqr).*abs(Ts_k).^2; 
HSK = fftshift(fft(HSK));
SARSpectrumPlot(abs(HSK),120,spectralBW,1);
title(['HSK Filter function']);

%% HintK
HintK = exp(-k_x.*xi_sqr).*(TR_k.*conj(Tvb_k) + conj(TR_k).*Tvb_k); 
HintK = fftshift(fft(HintK));
SARSpectrumPlot(abs(HintK),120,spectralBW,1);
title(['HintK Filter function']);
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
mean_VV = mean(sarData);
cvar = var((sarData-mean_VV)./mean_VV);
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
%%
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
%% Testing integrals
%dk = k(2) - k(1);
dth = theta(2) - theta(1);
int_th_trapz = cumtrapz(D(:,1))*dth;
figure;
plot(theta,int_th_trapz,'DisplayName' ,'Area under D(\theta)');
xlabel('\theta [rad]');
ylabel('\int D(\theta) \cdot d\theta')
grid on;
set(gca,'defaultAxesTickLabelInterpreter','latex'); 
set(gca,'XTick',L:pi/2:H) 
set(gca,'XTickLabel',{'','\theta_{wave}-\pi/2','\theta_{wave}','\theta_{wave}+\pi/2',''})
set(gca,'FontSize',12)
legend('show')
yline(1,'LineStyle','--')
matlab2tikz('../plots/directionalSpreadingVerify.tex');

