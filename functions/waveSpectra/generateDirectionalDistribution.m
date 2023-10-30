function [D,theta,s] = generateDirectionalDistribution(waveStruct,freqArray,freqChoice)

Tsig = waveStruct.significantWavePeriod(1,1);
Tp = Tsig./0.95;
switch freqChoice
    case 0
        freqSig = 1./Tsig;
        freqP = 1./Tp;
    case 1
        freqSig = 2*pi./Tsig;
        freqP = 2*pi./Tp;
end
sig_th = 26.9.*((freqSig)/(freqP)).^(-1.05); % in degrees
for i=1:length(freqArray)
    if (freqArray(i)>=freqP)
        sig_th(i) = 26.9.*((freqSig)/(freqP)).^(0.68);
    else
        sig_th(i) = 26.9.*((freqSig)/(freqP)).^(-1.05); % in degrees

    end
end
s = 2./(deg2rad(sig_th).^2) - 1;
Tp = waveStruct.significantWavePeriod(1,1)./0.95;
th_wave = waveStruct.direction(1,1);
th_wind = waveStruct.windDirection(1,1);

low = th_wave - pi;
high = th_wave + pi;
L = low  + rem(th_wave - low,0.025);
H = high - rem(high - th_wave,0.025);
theta = linspace(L,H,512)';

i = 1;
j = 1;
while (i<=512)
    A2 = gamma(s(i)+1)./(gamma(s(i)+1/2).*2.*sqrt(pi));
    D(:,j) = abs(A2 .* cos(0.5.*(th_wave-theta)).^(2*s(i))); % abs because D is complex. 
    %E(:,:,:,j) = S .* D(:,j)';
    i = i + 511;
    j = j + 1;
end


end