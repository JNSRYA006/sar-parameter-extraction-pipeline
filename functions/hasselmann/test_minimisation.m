%% Manual Method
func = helperFunctions;

%% Iterate cost function
numOfIterations = 7;

iterations = linspace(1,numOfIterations,numOfIterations);
sigWaveHeight = zeros(numOfIterations,1);
sigWavePeriod = zeros(numOfIterations,1);
% dk_x = 0.005890997229533;
% dk_y = 0.005226593472072;
dk_x = dkx;
dk_y = dky;
P_s_lin = func.resize(P_s_lin(:,2:end),E_k);
E_k = func.resize(E_k(2:end,2:end),P_s_lin);
E_k_inv = func.resize(E_k_inv(2:end,2:end),E_k);
Ts_k = func.resize(Ts_k(:,2:end),E_k);
Ts_k_inv = func.resize(Ts_k_inv(:,2:end),E_k);
p_s_coeff_inv = quasilinearCoeff(k_inv,k_y_inv,k_x_inv,E_k_inv,metadata,incidenceAngle);



for j=1:numOfIterations
numOfIterations = iterations(j);
[J(j),deltaEn,En] = costFunctionCalculation(P_s_pipeline,intensityFFT,E_k,E_k_inv,B,mu,Ts_k,Ts_k_inv,P_s_lin,p_s_coeff,p_s_coeff_inv,numOfIterations,dk_x,dk_y);
Pn = generateSARSpectrumOceanWaves(k,k_y,k_x,En,metadata,incidenceAngle,r,1,w);
En = abs(En);
%J_alt(j) = abs(costFunction2(P_s_pipeline,intensityFFT,E_k,En,B,mu,dk_x,dk_y));
deltaEn = func.resize(deltaEn(2:end,2:end),E_k);
En = func.resize(En(2:end,2:end),E_k);

% Convert back to th, w
d = 70;
c = sqrt((g./k).*tanh(k.*d)); % tanh in radians (eq. 5.4.23 in holthuisjen) output in m/s^2
n = 0.5*(1+(2.*k.*d)./sinh(2.*k.*d)); % sinh in radians (eq. 5.4.32 in holthuisjen)
c_g = n.*c;
En_w_th = En./((c.*c_g)./w);
En_w_th = func.resize(En_w_th(:,2:end),E_k);

for i=2:length(theta)
    dth(i) = theta(i)-theta(i-1);
end
dth = mean(dth);

%En_w = En_w_th./D(:,2);
En_w = trapz(En_w_th,2).*dth;

% Get out wave parameters
peakVal = max(En_w);
peakIndex = find(En_w == peakVal);
wPeak = w(peakIndex);
wSig = 0.95*wPeak;
sigWavePeriod(j) = (0.95*2*pi)./wPeak;
%sigWavePeriod(j) = (2*pi)./wPeak;

for i=2:length(w)
    dw(i) = w(i)-w(i-1);
end
dw = mean(dw);

int_En = trapz(trapz(En_w_th).*dth).*dw;
int_En_w = trapz(En_w).*dw;
sigWaveHeight(j) = 4.*sqrt(int_En_w);

disp(['Significant Wave Height = ', num2str(sigWaveHeight(j,1)), ' with iterations = ', num2str(numOfIterations)]);
disp(['Accuracy = ', num2str(sigWaveHeight(j,1)/waveVals.significantWaveHeight(2,2)), ' with iterations = ', num2str(numOfIterations)]);
disp(['Significant Wave Period = ', num2str(sigWavePeriod(j,1)), ' with iterations = ', num2str(numOfIterations)]);
disp(['Accuracy = ', num2str(sigWavePeriod(j,1)/waveVals.significantWavePeriod(2,2)), ' with iterations = ', num2str(numOfIterations)]);
end
%% Plot minimisation
figure;
plot(iterations,sigWaveHeight);
%% Plot output wave spectra

[t,r] = meshgrid(theta,w);
[x,y] = pol2cart(t,r);

figure;
%surf(x,y,En_w_th,'LineStyle','none');
contour(x,y,En_w_th);
title(['Best fit wave spectrum for number of iterations: ',num2str(numOfIterations)]);
yline(0);
xline(0);
grid on;
xlabel('\omega (rad/s)'), ylabel('\omega (rad/s)');
%%
EDiff = abs(En_w_th - E);
figure;
surf(x,y,EDiff,'LineStyle','none');
%contour(x,y,EDiff);
%title(['Number of iterations: ',num2str(numOfIterations)]);
yline(0);
xline(0);
grid on;
xlabel('\omega (rad/s)'), ylabel('\omega (rad/s)');
%%
figure;
%surf(k_x,k_y,En,'LineStyle','none')
contour(k_x,k_y,En);
title(['Best fit wave spectrum for number of iterations: ',num2str(numOfIterations)]);
grid on;
yline(0);
xline(0);
xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex');
%%
EDiff_k = abs(En - E_k);
figure;
%surf(k_x,k_y,EDiff_k,'LineStyle','none')
contour(k_x,k_y,EDiff_k);
% hold on;
% contour(k_x,k_y,E_k);
% contour(k_x,k_y,En);
% hold off;
%title(['Number of iterations: ',num2str(numOfIterations)]);
grid on;
yline(0);
xline(0);
xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex');
%% Plot output SAR spectra
figure;
SARSpectrumPlot(20*log10(abs(Pn)),135,spectralBW);
title('Minimised SAR Spectrum')

%%
figure;
%surf(k_x,k_y,20*log10(abs(Pn)),'LineStyle','none')
contour(k_x,k_y,20*log10(abs(Pn)));
title(['Best Fit SAR Spectrum for number of iterations: ',num2str(numOfIterations)]);
grid on;
yline(0);
xline(0);
xlabel('$k_{x}$','interpreter','latex'), ylabel('$k_{y}$','interpreter','latex');
%% Plot 1D spectrum
figure;
plot(w,En_w);
%% Plot 1D difference
E1DDiff = abs(En_w - S(:,:,1));
figure;
plot(w,E1DDiff,'DisplayName','Difference');
hold on;
plot(w,En_w,'DisplayName','Best-fit spectrum');
plot(w,S(:,:,1),'DisplayName','First-guess spectrum');
hold off;
legend('show');
%% Plot J
figure;
plot(iterations,J)


function [J,deltaFn,Fn] = costFunctionCalculation(generatedSARSpectrum,observedSARSpectrum,firstGuessWaveSpectrum,inverseFirstGuessWaveSpectrum,B,mu,Ts_k,Ts_k_inv,imageVarianceSpectrum,quasilinearCoeff,invQuasilinearCoeff,numOfIterations,dk_x,dk_y)

deltaFn = 0;
deltaPn = 0;
Fn = firstGuessWaveSpectrum;
Pn = observedSARSpectrum;
J = zeros(numOfIterations,1);

for i=1:numOfIterations
    Fn = Fn + deltaFn;
    Pn = Pn + deltaPn;
    deltaPk = observedSARSpectrum - Pn;
    deltaFk = firstGuessWaveSpectrum - Fn;
    deltaFnk = inverseFirstGuessWaveSpectrum - Fn;
    Wk = abs(Ts_k).^2.*quasilinearCoeff;
    Wnk = abs(Ts_k_inv).^2.*invQuasilinearCoeff;
    Ak = Wk.^2 + 2.*mu;
    Ank = Wnk.^2 + 2.*mu;
    Bk = Wk.*Wnk;
    deltaFn = (Ank.*(Wk.*deltaPk + mu.*deltaFk) - Bk.*(Wnk.*deltaPk + mu.*deltaFnk))./(Ak.*Ank - Bk.^2);
    deltaPn = 0.5.*quasilinearCoeff.*imageVarianceSpectrum;
    J = trapz(trapz((deltaPn-deltaPk).^2).*dk_x).*dk_y - mu.*trapz(trapz((deltaFn-deltaFk).^2).*dk_x).*dk_y;
end

end

function J = costFunction2(generatedSARSpectrum,observedSARSpectrum,firstGuessWaveSpectrum,optimalWaveSpectrum,B,mu,dk_x,dk_y)
% Equation 63 in HH
J = trapz(trapz((generatedSARSpectrum-observedSARSpectrum).^2.*observedSARSpectrum).*dk_x).*dk_y + mu.*trapz(trapz(((optimalWaveSpectrum - firstGuessWaveSpectrum)./(B + firstGuessWaveSpectrum)).^2).*dk_x).*dk_y;

end