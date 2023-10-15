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
% Update metadata function to general get look
%th = ncread(SARmetadata.Filename,"Incidence_Angle");

look = func.look(look);

% Need to check resizing is correct
%k_new = func.resize(k,th(1,:));
k_new = k;
%k_y = func.resize(k_y,th(1,:));
k_l = func.kl(look,k_y);
omega = func.omega(k_new);
%% Could be wrong
% Create logical masks for NaN rows and NaN columns
nanRows = any(isnan(waveSpectrum), 2); % For NaN rows
nanCols = any(isnan(waveSpectrum), 1); % For NaN columns

% Find row indices where NaN rows exist
nanRowIndices = find(nanRows);

% nanRowIndices will now contain the row indices of NaN rows
% Initialize masks of zeros with the same size as waveSpectrum
nanRowsMask = zeros(size(waveSpectrum));
nanColsMask = zeros(size(waveSpectrum));

% Set 1s in the masks where NaN values are in waveSpectrum
nanRowsMask(:,1) = 1;
nanColsMask(1,:) = 1;

% Now nanRowsMask and nanColsMask have 1s where NaN values are in waveSpectrum

% Remove NaN values manually
waveSpectrum = waveSpectrum(2:end,2:end);

% Resize the matrix to 2001x2001
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
% Figure out r value!!
%r = ones(2001);
%f_v_r = cumtrapz(waveSpectrum.*abs(Tv_k).^2).*dth;
f_v_r = cumtrapz(waveSpectrum.*abs(Tv_k).^2.*exp(1i.*k_new.*r)).*dk;

end



