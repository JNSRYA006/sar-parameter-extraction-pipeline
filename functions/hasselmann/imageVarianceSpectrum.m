function p_s_k = imageVarianceSpectrum(k,k_x,k_y,k_inv,k_x_inv,k_y_inv,waveSpectrum,waveSpectrum_inv,SARmetadata,th)
% Equation 26 in Hasselmann
%% Get required metadata
func = helperFunctions;
look = func.getLook(SARmetadata);
% Update metadata function to general get look
%th = ncread(SARmetadata.Filename,"Incidence_Angle");

look = func.look(look);
polarisation = func.getPolarisation(SARmetadata);
%captureTime = func.getCaptureTime(SARmetadata);
beta = func.getBeta(SARmetadata);

% Need to check resizing is correct
k_new = func.resize(k,th(1,:));
k_y = func.resize(k_y,th(1,:));
k_x = func.resize(k_x,th(1,:));
k_l = func.kl(look,k_y);
k_l_inv = func.kl(look,k_y_inv);
omega = func.omega(k_new);
mu=0.5;

Tt_k = func.tiltMTF(polarisation,k_l,th);
Th_k = func.hydroMTF(omega,mu,k_new,k_y);
Tt_k_inv = func.tiltMTF(polarisation,k_l_inv,th);
Th_k_inv = func.hydroMTF(omega,mu,k_inv,k_y_inv);
Tr_k = func.rarMTF(Tt_k,Th_k);
Tr_k_inv = func.rarMTF(Tt_k_inv,Th_k_inv);

Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tv_k_inv = func.rangeVelocityTF(omega,th,k_l_inv,k_inv);
Tvb_k = func.velocityBunchingMTF(beta,k_x,Tv_k);
Tvb_k_inv = func.velocityBunchingMTF(beta,k_x_inv,Tv_k_inv);

Ts_k = func.sarImagingMTF(Tr_k,Tvb_k);
Ts_k_inv = func.sarImagingMTF(Tr_k_inv,Tvb_k_inv);

% Remove NaN values manually
waveSpectrum = waveSpectrum(2:end,2:end);

% Resize the matrix to 2001x2001
waveSpectrum = func.resize(waveSpectrum,th);

% Remove NaN values manually
waveSpectrum_inv = waveSpectrum_inv(2:end,1:end-1);

% Resize the matrix to 2001x2001
waveSpectrum_inv = func.resize(waveSpectrum_inv,th);
    
p_s_k = abs(Ts_k).^2.*(waveSpectrum./2) + abs(Ts_k_inv).^2.*(waveSpectrum_inv./2);
end