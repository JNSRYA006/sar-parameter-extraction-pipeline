function B = calculateNormalisationConstant(firstGuessWaveSpectrum)
B = 0.01.*max(firstGuessWaveSpectrum(:));
end