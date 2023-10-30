function [f_r_r] = rarImageIntensityAutocovariance(k,k_y,waveSpectrum,SARmetadata,r,th)
    
% Need:
% - F(k)
% - T^R_k
% -- T^t_k
%   - k_l (Wavenumber in look direction)
%   -- Look
%   - theta (Radar incidence angle)
%   - Polarisation
% -- T^h_k
%   - omega
%   - mu
%   - k
%   - k_y
% -- theta (Radar incidence angle)
% -- k_l (Wavenumber in look direction)
% - k and r (vectors)
% -- r is the x and y coords
% function: f = @(x,y) (x^.2 + y.^2)

% x-axis = SAR flight direction (azimuthal) (heading)

func = helperFunctions;

look = func.getLook(SARmetadata);
polarisation = func.getPolarisation(SARmetadata);

look = func.look(look);
k_new = k;
k_l_inv = func.kl(look,fliplr(-k_y));
k_l = func.kl(look,k_y);
omega = func.omega(k_new);

waveSpectrum = func.resize(waveSpectrum,th);
mu = 0.5;

Tt_k = func.tiltMTF(polarisation,k_l,th);
Th_k = func.hydroMTF(omega,mu,k_new,k_y);
Th_k = func.resize(Th_k(2:end),waveSpectrum(1,:));
Tt_k_inv = func.tiltMTF(polarisation,k_l_inv,th);
Th_k_inv = func.hydroMTF(omega,mu,fliplr(-k_new),k_y);
Th_k_inv = func.resize(Th_k_inv(1:end-1),waveSpectrum(1,:));
Tr_k = func.rarMTF(Tt_k,Th_k);
Tr_k_inv = func.rarMTF(Tt_k_inv,Th_k_inv);

for i=2:length(k)
    dk(i) = k(i)-k(i-1);
end
dk = dk(2:end);
dk = mean(dk);

f_r_r = 0.5.*cumtrapz((waveSpectrum.*abs(Tr_k).^2)+(waveSpectrum.*abs(Tr_k_inv).^2).*exp(1i.*k_new.*r)).*dk;


end