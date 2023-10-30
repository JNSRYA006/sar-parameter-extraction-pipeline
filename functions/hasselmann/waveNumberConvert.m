function [E_w, E_w_th] = waveNumberConvert(En,k,w,theta,D,d)
    % Convert back to th, w
    %d = 70;
    g = 9.81;
    c = sqrt((g./k).*tanh(k.*d)); % tanh in radians (eq. 5.4.23 in holthuisjen) output in m/s^2
    n = 0.5*(1+(2.*k.*d)./sinh(2.*k.*d)); % sinh in radians (eq. 5.4.32 in holthuisjen)
    c_g = n.*c;
    En_w_th = En./((c.*c_g)./w);
    En_w_th = func.resize(En_w_th(:,2:end),E_k);
    
    for i=2:length(theta)
        dth(i) = theta(i)-theta(i-1);
    end
    dth = mean(dth);
    
    En_w = En_w_th/D(:,1)';
    %En_w = trapz(En_w_th,2).*dth;
end