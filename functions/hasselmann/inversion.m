function [J,En,Pn,En_w_th,En_w,sigWaveHeight,sigWavePeriod] = inversion(iterations,generatedSARSpectrum,observedSARSpectrum,firstGuessWaveSpectrum,inverseFirstGuessWaveSpectrum,k,k_x,k_y,B,mu,Ts_k,Ts_k_inv,k_x_inv,k_y_inv,imageVarianceSpectrum,quasilinearCoeff,invQuasilinearCoeff,dk_x,dk_y,metadata,incidenceAngle,r,w,D,theta)

%% Iterate cost function
func = helperFunctions;
g = 9.81;
iterationsLin = linspace(1,iterations,iterations);
sigWaveHeight = zeros(iterations,1);
sigWavePeriod = zeros(iterations,1);
% dk_x = 0.005890997229533;
% dk_y = 0.005226593472072;
P_s_lin = func.resize(imageVarianceSpectrum(:,2:end-1),firstGuessWaveSpectrum);
E_k = func.resize(firstGuessWaveSpectrum(2:end,2:end),P_s_lin);
E_k_inv = func.resize(inverseFirstGuessWaveSpectrum(2:end,1:end-1),firstGuessWaveSpectrum);
Ts_k = func.resize(Ts_k(:,2:end),firstGuessWaveSpectrum);
Ts_k_inv = func.resize(Ts_k_inv(:,1:end-1),firstGuessWaveSpectrum);

% p_s_coeff_inv = quasilinearCoeff(k_inv,k_y_inv,k_x_inv,inverseFirstGuessWaveSpectrum,metadata,incidenceAngle);
% 
% 
Pn = generatedSARSpectrum;
En = firstGuessWaveSpectrum;

for j=1:iterations
    numOfIterations = iterationsLin(j);
    [J(j),deltaEn,En,deltaPn] = costFunctionCalculation(generatedSARSpectrum,observedSARSpectrum,En,E_k_inv,B,mu,Ts_k,Ts_k_inv,P_s_lin,quasilinearCoeff,invQuasilinearCoeff,numOfIterations,dk_x,dk_y);
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
    
    En_w = En_w_th/D(:,1)';
    %En_w = trapz(En_w_th,2).*dth;
    
    % Get out wave parameters
    peakVal = max(En_w);
    peakIndex = find(En_w == peakVal);
    wPeak = w(peakIndex);
    wSig = 0.95*wPeak;
    %sigWavePeriod(j) = (2*pi)./wSig;
    sigWavePeriod(j) = (0.95*2*pi)./wPeak;
    %sigWavePeriod(j) = (2*pi)./wPeak;
    
    for i=2:length(w)
        dw(i) = w(i)-w(i-1);
    end
    dw = mean(dw);
    
    int_En = trapz(trapz(En_w_th).*dth).*dw;
    int_En_w = trapz(En_w).*dw;
    sigWaveHeight(j) = 4.*sqrt(int_En_w);
    % heightAcc(j) = sigWaveHeight(j,1)/waveVals.significantWaveHeight(2,2);
    % periodAcc(j) = sigWavePeriod(j,1)/waveVals.significantWavePeriod(2,2);
    % 
    disp(['Significant Wave Height = ', num2str(sigWaveHeight(j,1)), ' with iterations = ', num2str(numOfIterations)]);
    % disp(['Accuracy = ', num2str(sigWaveHeight(j,1)/waveVals.significantWaveHeight(2,2)), ' with iterations = ', num2str(numOfIterations)]);
    disp(['Significant Wave Period = ', num2str(sigWavePeriod(j,1)), ' with iterations = ', num2str(numOfIterations)]);
    % disp(['Accuracy = ', num2str(sigWavePeriod(j,1)/waveVals.significantWavePeriod(2,2)), ' with iterations = ', num2str(numOfIterations)]);
end
end


function J = costFunction(generatedSARSpectrum,observedSARSpectrum,firstGuessWaveSpectrum,B,mu,dk)
% Equation 63 in HH
J = cumtrapz((generatedSARSpectrum-observedSARSpectrum).^2).*observedSARSpectrum.*dk + mu.*cumtrapz(((optimalWaveSpectrum - firstGuessWaveSpectrum)./(B + firstGuessWaveSpectrum)).^2).*dk;

end