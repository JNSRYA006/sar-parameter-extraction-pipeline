function [f_r_r] = rarImageIntensityAutocovariance(k,k_y,waveSpectrum,SARmetadata,th,r)
    
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

%% Get required metadata
[look,inc_near,inc_far,num_pixels,polarisation] = getMetadata_f_r_r(SARmetadata);

func = helperFunctions;

look = func.look(look);



%th = func.incidence(inc_near,inc_far,num_pixels);
%th = extrapolateIncidence(inc_near,inc_far,num_pixels); % Need to calculate based on velocity and look time

% Need to check resizing is correct
k_new = func.resize(k,th(1,:));
k_y = func.resize(k_y,th(1,:));
k_l = func.kl(look,k_y);
omega = func.omega(k_new);
%% Could be wrong
% Create a logical mask for rows with NaN entries
nanRows = any(isnan(waveSpectrum), 2);

% Remove rows with NaN entries
waveSpectrum = waveSpectrum(~nanRows, :);

waveSpectrum = func.resize(waveSpectrum,zeros(2001,2001));
mu = 0;

Tt_k = func.tiltMTF(polarisation,k_l,th);
Th_k = func.hydroMTF(omega,mu,k_new,k_y);
Tt_k_inv = func.tiltMTF(polarisation,k_l,th);
Th_k_inv = func.hydroMTF(omega,mu,-k_new,k_y);
Tr_k = Tt_k + Th_k;
Tr_k_inv = Tt_k_inv + Th_k_inv;

% k_new = resizeToSameSize(k,th);
% k_y = resizeToSameSize(k_y,th);
% k_l = defineKLook(look,k_y);
% omega = defineOmega(k_new);
% waveSpectrum = resizeToSameSize(waveSpectrum,zeros(2001,2001));

for i=2:length(k)
    dk(i) = k(i)-k(i-1);
end
dk = dk(2:end);
dk = mean(dk);
% Figure out r value!!
%r = ones(2001);

% Need to figure out how to get F(-k)

f_r_r = 0.5.*cumtrapz((waveSpectrum.*abs(Tr_k).^2)+(waveSpectrum.*abs(Tr_k_inv).^2).*exp(1i.*k_new.*r)).*dk;


end

function [look,incidence_near,incidence_far,num_pixels,polarisation] = getMetadata_f_r_r(metadata)
    req_atributes = ["antenna_pointing","incidence_near","incidence_far","num_output_lines","mds1_tx_rx_polar","mds2_tx_rx_polar"];
    meta_frr = filterAttributesNetCDF(metadata.Attributes, req_atributes);
    polarisation= ["";""];
    look = meta_frr(1).Value;
    incidence_near = meta_frr(2).Value;
    incidence_far = meta_frr(3).Value;
    num_pixels = meta_frr(4).Value;
    polarisation(1) = meta_frr(5).Value;
    polarisation(2) = meta_frr(6).Value;
end
