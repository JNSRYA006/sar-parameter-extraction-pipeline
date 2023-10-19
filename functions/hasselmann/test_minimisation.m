%% Manual Method
func = helperFunctions;

%% Iterate cost function
numOfIterations = 6;

iterations = linspace(1,numOfIterations,numOfIterations);
sigWaveHeight = zeros(numOfIterations,1);
sigWavePeriod = zeros(numOfIterations,1);
dk_x = 0.005890997229533;
dk_y = 0.005226593472072;
P_s_lin = func.resize(P_s_lin(:,2:end),E_k);
E_k = func.resize(E_k(2:end,2:end),P_s_lin);
E_k_inv = func.resize(E_k_inv(2:end,2:end),E_k);
Ts_k = func.resize(Ts_k(:,2:end),E_k);
Ts_k_inv = func.resize(Ts_k_inv(:,2:end),E_k);
p_s_coeff_inv = quasilinearCoeff(k_inv,k_y_inv,k_x_inv,E_k_inv,metadata,incidenceAngle);

for j=1:numOfIterations
numOfIterations = iterations(j);
%numOfIterations = 8;
[J,deltaEn,En] = costFunctionCalculation(P_s_pipeline,intensityFFT,E_k,E_k_inv,B,mu,Ts_k,Ts_k_inv,P_s_lin,p_s_coeff,p_s_coeff_inv,numOfIterations,dk_x,dk_y);
Pn = generateSARSpectrumOceanWaves(k,k_y,k_x,En,metadata,incidenceAngle,r,1,w);
En = abs(En);

deltaEn = func.resize(deltaEn(2:end,2:end),E_k);
En = func.resize(En(2:end,2:end),E_k);

% Convert back to th, w
d = 70;
c = sqrt((g./k).*tanh(k.*d)); % tanh in radians (eq. 5.4.23 in holthuisjen) output in m/s^2
n = 0.5*(1+(2.*k.*d)./sinh(2.*k.*d)); % sinh in radians (eq. 5.4.32 in holthuisjen)
c_g = n.*c;
En_w_th = En./((c.*c_g)./w);
En_w_th = func.resize(En_w_th(:,2:end),E_k);

%En_w = En_w_th./D(:,2);
En_w = trapz(En_w_th,2).*dth;

% Get out wave parameters
peakVal = max(En_w);
peakIndex = find(En_w == peakVal);
wPeak = w(peakIndex);
wSig = 0.95*wPeak;
%sigWavePeriod(j) = (0.95*2*pi)./wPeak;
sigWavePeriod(j) = (2*pi)./wPeak;

for i=2:length(theta)
    dth(i) = theta(i)-theta(i-1);
end
dth = mean(dth);

for i=2:length(w)
    dw(i) = w(i)-w(i-1);
end
dw = mean(dw);

int_En = trapz(trapz(En).*dth).*dw;
int_En_w = trapz(En_w).*dw;
sigWaveHeight(j) = 4.*sqrt(int_En);

disp(['Significant Wave Height = ', num2str(sigWaveHeight(j,1)), ' with iterations = ', num2str(numOfIterations)]);
disp(['Significant Wave Period = ', num2str(sigWavePeriod(j,1)), ' with iterations = ', num2str(numOfIterations)]);
end
%% Plot minimisation
figure;
plot(iterations,sigWavePeriod);
%% Plot output wave spectra

[t,r] = meshgrid(theta,w);
[x,y] = pol2cart(t,r);

figure;
%surf(x,y,En_w_th,'LineStyle','none');
contour(x,y,En_w_th);
%title(['Number of iterations: ',num2str(numOfIterations)]);
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
%title(['Number of iterations: ',num2str(numOfIterations)]);
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
surf(k_x,k_y,20*log10(abs(Pn)),'LineStyle','none')
%contour(k_x,k_y,20*log10(abs(Pn)));
%title(['Number of iterations: ',num2str(numOfIterations)]);
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
plot(w,E1DDiff);
hold on;
plot(w,En_w);
plot(w,S(:,:,1))
hold off;
%% Plot J
figure;
plot(iterations,abs(J))


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
    J(i) = trapz(trapz((deltaPn-deltaPk).^2).*dk_x).*dk_y + mu.*trapz(trapz((deltaFn-deltaFk).^2).*dk_x).*dk_y;
end

end

function J = costFunction(generatedSARSpectrum,observedSARSpectrum,firstGuessWaveSpectrum,B,mu,dk)
% Equation 63 in HH
J = cumtrapz((generatedSARSpectrum-observedSARSpectrum).^2).*observedSARSpectrum.*dk + mu.*cumtrapz(((optimalWaveSpectrum - firstGuessWaveSpectrum)./(B + firstGuessWaveSpectrum)).^2).*dk;

end
%% Tried using fmincon
% dk = 0.015991335372069;
% %optimalWaveSpectrum = zeros(2001,2001);
% % Define your cost function as a function of optimalWaveSpectrum
% costFunction = @(optimalWaveSpectrum) computeCost(optimalWaveSpectrum,abs(P_s),VV_nc,E_k_resize,B,mu,dk_x,dk_y);
% 
% % Initial guess for optimalWaveSpectrum
% func = helperFunctions;
% % Remove rows with NaN entries
% E_k_resize = E_k(2:end,2:end);
% E_k_resize = func.resize(E_k_resize,th);
% %initialGuess = E_k;
% 
% %test = computeCost(optimalWaveSpectrum,abs(P_s),VV_nc,E_k,B,mu,dk_x,dk_y);
% 
% % Define initial guess
% initialGuess = zeros(512, 512);  % Initialize with zeros or other appropriate values
% 
% % Define lower and upper bounds for each element of optimalWaveSpectrum
% lb = -Inf(512, 512);  % Lower bound (no constraint)
% ub = Inf(512, 512);   % Upper bound (no constraint)
% 
% % Define nonlinear constraints (empty in this case)
% nonlcon = [];
% 
% % Define options for fmincon
% options = optimoptions(@fmincon, 'Display', 'iter', 'OutputFcn', @outputFunction);  % Add custom output function
% 
% % Perform constrained optimization
% [optimalWaveSpectrum, minimumCost] = fmincon(costFunction, initialGuess, [], [], [], [], lb, ub, nonlcon, options);
% 
% % Display the optimized result
% disp(['Optimal Wave Spectrum: ', num2str(optimalWaveSpectrum)]);
% disp(['Minimum Cost: ', num2str(minimumCost)]);
% 
% % Define the cost function
% function cost = computeCost(optimalWaveSpectrum, generatedSARSpectrum, observedSARSpectrum, firstGuessWaveSpectrum, B, mu, dk_x, dk_y)
%     % Calculate the difference between generatedSARSpectrum and observedSARSpectrum
% diff1 = generatedSARSpectrum - observedSARSpectrum;
% integral1 = trapz(trapz(diff1.^2).*dk_x).* dk_y;
% 
% % Calculate the difference between optimalWaveSpectrum and firstGuessWaveSpectrum
% diff2 = optimalWaveSpectrum - firstGuessWaveSpectrum;
% integral2 = trapz(trapz((diff2./(B + firstGuessWaveSpectrum)).^2).*dk_x).* dk_y;
% 
% % Combine the integrals with appropriate weights
% cost = integral1 + mu * integral2;
% 
% end
% 
% % Custom output function to display the number of iterations
% function stop = outputFunction(x, optimValues, state)
%     stop = false;
%     if isequal(state, 'iter')
%         disp(['Iteration: ' num2str(optimValues.iteration)]);
%     end
% end
