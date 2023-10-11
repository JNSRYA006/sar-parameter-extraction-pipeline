function E = generate2DWaveSpectrum(waveSpectrum,D)
i = 1;
j = 1;
while (i<=252)
    E(:,:,:,j) = waveSpectrum .* D(:,j)';
    i = i + 251;
    j = j + 1;
end

end