function SARSpectrumPlot(intensityFFT,spectrumThreshold,spectralBW,normaliseAxes)

    spectrumIndex = intensityFFT > spectrumThreshold;
    intensityFFT = spectrumIndex.*intensityFFT;
    
    % Calculate the spatial frequencies
    [M, N] = size(intensityFFT);
    dx = spectralBW;
    dy = spectralBW;
    kx = (-M/2:M/2-1) / (M * dx);
    ky = (-N/2:N/2-1) / (N * dy);
    
    % Calculate Δk (the spacing between k values)
    dkx = kx(2) - kx(1);
    dky = ky(2) - ky(1);
    
    subIntensityFFT = intensityFFT(231:281,231:281);
    % subIntensityFFT = intensityFFT;
    if normaliseAxes
        % Create a grid for the contour plot
        [X, Y] = meshgrid(kx/dkx, ky/dky);
        X = X(231:281,231:281);
        Y = Y(231:281,231:281);
        labelX = 'k_x/\Deltak';
        labelY = 'k_y/\Deltak';
    else
        [X, Y] = meshgrid(kx, ky);
        X = X(231:281,231:281);
        Y = Y(231:281,231:281);
        labelX = 'k_x';
        labelY = 'k_y';
    end
    % Plot the SAR spectrum with k_x/Δk and k_y/Δk axes using contourf
    figure;
    %contourf(X, Y, 20*log10(abs(subIntensityFFT)), 20, 'LineColor', 'none');
    %surf(X,Y,20*log10(subIntensityFFT),'LineStyle','none');
    contour(X, Y, 20*log10(subIntensityFFT),'LineWidth',0.5);
    %imagesc(kx/dkx,ky/dky,20*log10(intensityFFT));
    colormap('jet');
    cb = colorbar;  % create and label the colorbar
    cb.Label.String = 'Power Spectrum [dB]';
    yline(0,'LineWidth',0.25);
    xline(0);
    grid on;
    grid minor;
    
    % % Set the X and Y axis ticks and labels
    % xticks(unique([min(kxSubset):1:max(kxSubset)]));
    % yticks(unique([min(kySubset):1:max(kySubset)]));
    % xticklabels(unique(string(kxSubset)));
    % yticklabels(unique(string(kySubset)));
    
    xlabel(labelX);
    ylabel(labelY);
    
    %title('Observed SAR Spectrum');


end