function E = generateMultipleJONSWAP(waveStruct,gammaVal,freqArray,freqChoice)

lat = waveStruct.latitude;
lon = waveStruct.longitude;

count = 1;
for i = 1:length(lat)
    for j = 1:length(lon)
        Hs = waveStruct.significantWaveHeight(j,i);
        T0 = waveStruct.significantWavePeriod(j,i);
        switch freqChoice
            case 0
                freq0 = 1./T0;
            case 1
                freq0 = 2*pi./T0;
        end
        E(:,:,count) = wavespec(7,[Hs,freq0,gammaVal],freqArray,0);
        count = count + 1;
    end
end

end