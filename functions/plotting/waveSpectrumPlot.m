function waveSpectrumPlot(waveVals,S,w,singleOrMultiple,displayw,plotAllLats,plotAllLons)

if singleOrMultiple
    figure;
    %plot(w,S/(T0*Hs^2))
    %title('S(\omega)/(H_s^2 T_0)')
    
    hold on;
    if plotAllLats
        n = length(waveVals.longitude);        
    else
        n=1;
    end
    if plotAllLons
        m = length(waveVals.latitude);
    else
        m=1;
    end
    for i = 1:n*m
        lat_index = ceil(i/m);
        lon_index = mod(i-1,m)+1;
        display_name = [num2str(waveVals.latitude(lat_index)), 'S, ' num2str(waveVals.longitude(lon_index)), 'E'];
        % For plotting multiple lat vals
        switch lat_index
            case 1
                colourToPlot = "#0072BD";
            case 2
                colourToPlot = "#EDB120";
            case 3 
                colourToPlot = "#A2142F";
        end
        plot(w,S(:,i),'DisplayName',display_name, 'Color',colourToPlot);
        %plot(w,S(:,i),'DisplayName',display_name);
        maxIndex = find(S(:,i) == max(S(:,i)));
        switch lon_index
            case 1
                colourToPlot = "#0072BD";
            case 2              
                colourToPlot = "#D95319";
        end
        S(maxIndex,i);
        if displayw
            xline(w(maxIndex),LineWidth=1,DisplayName=['\omega = ',num2str(round(w(maxIndex),3))],Color=colourToPlot,LineStyle="--");
        end
        %legend(["Wave spectrum at: ", num2str(testLat(lat_index)), ", " num2str(testLon(lon_index))]);
        %fprintf('Iteration %d: Latitude %.2f, Longitude %.2f\n', i, testLat(lat_index), testLon(lon_index));
    end
    hold off

    %figure;
    %plot(w,S(:,1));
    %plot(f,S);
    title('One-dimensional wave spectrum, E(\omega)')
    %legend('show');
    %legend('JONSWAP for \gamma = 1.308','Location','Northeast')
    xlabel('\omega [rad/s]')
    ylabel('E(\omega) [m^2/rad/Hz]')
    grid on;
    % set(gca,'XTick',0:0.25:1) 
    % set(gca,'XTickLabel',{'0','0.25','0.5','0.75','1'})
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'FontSize',12)
    legend('show')
    %matlab2tikz('../plots/valWaveSpec_toShore_1.tex');
else
    figure;
    plot(w,S);
    %plot(f,S);
    title(['One-dimensional wave spectrum, E(\omega) at ', num2str(waveVals.latitude(1,1)), 'S, ', num2str(waveVals.longitude(1,1)), 'E'])
    legend('JONSWAP for \gamma = 1.308','Location','Northeast')
    xlabel('\omega [rad/s]')
    ylabel('E(\omega) [m^2/rad/Hz]')
    grid on;
    % set(gca,'XTick',0:0.25:1) 
    % set(gca,'XTickLabel',{'0','0.25','0.5','0.75','1'})
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'FontSize',12)
    legend('show')
    %matlab2tikz('../plots/oneDWaveSpec.tex');
end    


end