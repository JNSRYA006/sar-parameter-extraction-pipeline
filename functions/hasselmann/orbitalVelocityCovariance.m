function f_v_r = orbitalVelocityCovariance(waveSpectrum,radarIncidenceAngle,look,heading)

% Need:
% - F(k)
% - T^V_k
% -- omega
% -- theta (Radar incidence angle)
% -- k_l (Wavenumber in look direction
% --- Look
% -- k
% - k and r (vectors)
% -- r is the x and y coords
% function: f = @(x,y) (x^.2 + y.^2)

% x-axis = SAR flight direction (azimuthal) (heading)

%% Define functions used
omega = @(k) (sqrt(g.*k));
Tv_k = @(omega,th,k_l,k) (-omega.*(sin(th)*k_l./abs(k) + 1i.*cos(th)));


% %% Single heading value (in degrees)
% heading = deg2rad(90 - 344.3488); % Replace with your desired heading angle
% 
% % Create a unit vector based on the heading angle
% u = cos(heading); % Calculate x-component of the unit vector
% v = sin(heading); % Calculate y-component of the unit vector
% quiver(17.6424, -34.0603, u, v, 'bl', 'LineWidth', 2); % Assuming origin at (centre_lon, centre_lat)

f_v_r = 0;
end