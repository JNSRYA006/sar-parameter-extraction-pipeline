function [E_k,k] = waveNumberSpectrum(waveSpectrum,w,k,d)

g = 9.81;
c = sqrt((g./k).*tanh(k.*d)); % tanh in radians (eq. 5.4.23 in holthuisjen) output in m/s^2
n = 0.5*(1+(2.*k.*d)./sinh(2.*k.*d)); % sinh in radians (eq. 5.4.32 in holthuisjen)
c_g = n.*c;
E_k = ((c.*c_g)./w).*waveSpectrum;

end