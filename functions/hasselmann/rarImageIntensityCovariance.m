function [f_rv_r] = rarImageIntensityCovariance(k,k_y,waveSpectrum,SARmetadata,r,th)
    
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
% -- T^V_k
%   - Same as above
% - k and r (vectors)
% -- r is the x and y coords
% function: f = @(x,y) (x^.2 + y.^2)

% x-axis = SAR flight direction (azimuthal) (heading)

func = helperFunctions;

look = func.getLook(SARmetadata);
polarisation = func.getPolarisation(SARmetadata);

look = func.look(look);
k_new = k;
k_l = func.kl(look,k_y);
k_l_inv = func.kl(look,fliplr(-k_y));
omega = func.omega(k_new);

waveSpectrum = waveSpectrum(2:end,2:end);

waveSpectrum = func.resize(waveSpectrum,th);
mu = 0.5;

Tt_k = func.tiltMTF(polarisation,k_l,th);
Th_k = func.hydroMTF(omega,mu,k_new,k_y);
Th_k = func.resize(Th_k(2:end),waveSpectrum(1,:));
Tt_k_inv = func.tiltMTF(polarisation,k_l_inv,th);
Th_k_inv = func.hydroMTF(omega,mu,fliplr(-k_new),fliplr(-k_y));
Th_k_inv = func.resize(Th_k_inv(1:end-1),waveSpectrum(1,:));
Tr_k = func.rarMTF(Tt_k,Th_k);
Tr_k_inv = func.rarMTF(Tt_k_inv,Th_k_inv);

Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tv_k = func.resize(Tv_k(:,2:end),waveSpectrum);
Tv_k_inv = func.rangeVelocityTF(omega,th,k_l_inv,fliplr(-k_new));
Tv_k_inv = func.resize(Tv_k_inv(:,1:end-1),waveSpectrum);

for i=2:length(k)
    dk(i) = k(i)-k(i-1);
end
dk = dk(2:end);
dk = mean(dk);

f_rv_r = 0.5.*cumtrapz((waveSpectrum.*Tr_k.*conj(Tv_k))+(waveSpectrum.*conj(Tr_k_inv).*Tv_k_inv).*exp(1i.*k_new.*r)).*dk;


end

function [look,incidence_near,incidence_far,num_pixels,polarisation] = getMetadata_f_rv_r(metadata)
    req_atributes = ["antenna_pointing","incidence_near","incidence_far","num_output_lines","mds1_tx_rx_polar","mds2_tx_rx_polar"];
    meta_frvr = filterAttributesNetCDF(metadata.Attributes, req_atributes);
    polarisation= ["";""];
    look = meta_frvr(1).Value;
    incidence_near = meta_frvr(2).Value;
    incidence_far = meta_frvr(3).Value;
    num_pixels = meta_frvr(4).Value;
    polarisation(1) = meta_frvr(5).Value;
    polarisation(2) = meta_frvr(6).Value;
end
