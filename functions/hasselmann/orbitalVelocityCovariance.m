function [f_v_r] = orbitalVelocityCovariance(k,k_y,waveSpectrum,SARmetadata,th,r)

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
[look,inc_near,inc_far,num_pixels] = getMetadata_f_v_r(SARmetadata);
% Update metadata function to general get look
%th = ncread(metadata.Filename,"Incidence_Angle"); % Check this works when
%I have hardrive
func = helperFunctions;

%look = lookDiscretise(look);
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
Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);

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

f_v_r = cumtrapz(waveSpectrum.*abs(Tv_k).^2.*exp(1i.*k_new.*r)).*dk;

%f_v_r = integral(waveSpectrum*abs(Tv_k(omega,th,k_l,k)).^2.*exp(1i.*k.*r));
% %% Single heading value (in degrees)
% heading = deg2rad(90 - 344.3488); % Replace with your desired heading angle
% 
% % Create a unit vector based on the heading angle
% u = cos(heading); % Calculate x-component of the unit vector
% v = sin(heading); % Calculate y-component of the unit vector
% quiver(17.6424, -34.0603, u, v, 'bl', 'LineWidth', 2); % Assuming origin at (centre_lon, centre_lat)

%f_v_r = 0;
end

function [look,incidence_near,incidence_far,num_pixels] = getMetadata_f_v_r(metadata)
    req_atributes = ["antenna_pointing","incidence_near","incidence_far","num_output_lines"];
    meta_fvr = filterAttributesNetCDF(metadata.Attributes, req_atributes);
    look = meta_fvr(1).Value;
    incidence_near = meta_fvr(2).Value;
    incidence_far = meta_fvr(3).Value;
    num_pixels = meta_fvr(4).Value;
end


