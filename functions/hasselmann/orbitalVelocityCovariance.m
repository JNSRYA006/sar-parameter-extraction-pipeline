function [f_v_r] = orbitalVelocityCovariance(k,k_y,waveSpectrum,SARmetadata,r,th)

% Need:
% - F(k)
% - T^V_k
% -- omega
% -- theta (Radar incidence angle)
% -- k_l (Wavenumber in look direction)
% --- Look
% -- k
% - k and r (vectors)
% -- r is the x and y coords
% function: f = @(x,y) (x^.2 + y.^2)

% x-axis = SAR flight direction (azimuthal) (heading)

%% Get required metadata
func = helperFunctions;
look = func.getLook(SARmetadata);

look = func.look(look);

k_new = k;
k_l = func.kl(look,k_y);
omega = func.omega(k_new);

% Remove NaN values manually
waveSpectrum = waveSpectrum(2:end,2:end);
waveSpectrum = func.resize(waveSpectrum,th);

Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tv_k = func.resize(Tv_k(:,2:end),waveSpectrum);

for i=2:length(k)
    dk(i) = k(i)-k(i-1);
end
dk = dk(2:end);
dk = mean(dk);

for i=2:length(th)
    dth(i) = th(i)-th(i-1);
end
dth = dth(2:end);
dth = mean(dth);

f_v_r = cumtrapz(waveSpectrum.*abs(Tv_k).^2.*exp(1i.*k_new.*r)).*dk;

end



