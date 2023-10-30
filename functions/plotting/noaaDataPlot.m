function noaaDataPlot(projection,waveStruct,dataTypeToPlot)
    

    lat = waveStruct.latitude;
    lon = waveStruct.longitude;
    switch dataTypeToPlot
        case 'significantWaveHeight'
            titleStr = ['Significant Wave Height from NOAA (NCEP) - ',char(waveStruct.time)];
            data = waveStruct.significantWaveHeight;
            barStr = '[m]';
        case 'significantWavePeriod'
            titleStr = ['Significant Wave Period from NOAA (NCEP) - ',char(waveStruct.time)];
            data = waveStruct.significantWavePeriod;
            barStr = '[s]';
        case 'direction'
            titleStr = ['Wave Direction from NOAA (NCEP) - ',char(waveStruct.time)];
            data = waveStruct.direction;
            barStr = '[degrees]';
        case 'windSpeed'
            titleStr = ['Wind Speed from NOAA (NCEP) - ',char(waveStruct.time)];
            data = waveStruct.windSpeed;
            barStr = '[m/s]';
        case 'windDirection'
            titleStr = ['Wind Direction from NOAA (NCEP) - ',char(waveStruct.time)];
            data = waveStruct.windDirection;
            barStr = '[degrees]';
    end

    fig = figure;
    m_proj(projection,'lat',[min(lat(:)) max(lat(:))],...
    'lon',[min(lon(:)) max(lon(:))])
    % Next, plot the field using the M_MAP version of pcolor.
    m_pcolor(lon,lat,data.');
    shading flat;
    % Add a coastline and axis values.
    m_coast('patch',[.7 .7 .7])
    m_grid('box','fancy')
    % Add a colorbar and title.
    c = colorbar;
    %c.Label.String = barStr;
    hL = ylabel(c,barStr);     
    set(hL,'Rotation',0);
    %orient(fig,'landscape')
    %print(fig,['../plots/noaa',projection,dataTypeToPlot,'.pdf'],'-dpdf','-fillpage')
    title(titleStr);
    %matlab2tikz(['../plots/noaa',projection,dataTypeToPlot,'.tex']);

end