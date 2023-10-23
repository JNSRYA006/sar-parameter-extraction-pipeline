function J = costFunction(generatedSARSpectrum,observedSARSpectrum,firstGuessWaveSpectrum,B,mu,dk)
% Equation 63 in HH
J = cumtrapz((generatedSARSpectrum-observedSARSpectrum).^2).*observedSARSpectrum.*dk + mu.*cumtrapz(((optimalWaveSpectrum - firstGuessWaveSpectrum)./(B + firstGuessWaveSpectrum)).^2).*dk;

end