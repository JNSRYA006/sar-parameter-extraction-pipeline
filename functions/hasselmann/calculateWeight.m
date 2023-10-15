function mu0 = calculateWeight(observedSARSpectrum)
mu0 = 0.1.*max(observedSARSpectrum(:)).^2;
end