function p_s_ql = quasilinearApprox(k,k_y,k_x,waveSpectrum,SARmetadata,th)

func = helperFunctions;

[beta,look] = getBeta(SARmetadata,245,17,'max');

look = func.look(look);
%th = func.incidence(inc_near,inc_far,num_pixels);
%th = extrapolateIncidence(inc_near,inc_far,num_pixels); % Need to calculate based on velocity and look time
% Need to check resizing is correct
k_new = func.resize(k,th);
k_y = func.resize(k_y,th);
k_x = func.resize(k_x,th);
k_l = func.kl(look,k_y);
omega = func.omega(k_new);
%% Could be wrong
% Create a logical mask for rows with NaN entries
nanRows = any(isnan(waveSpectrum), 2);

% Remove rows with NaN entries
waveSpectrum = waveSpectrum(~nanRows, :);

waveSpectrum = func.resize(waveSpectrum,zeros(2001,2001));
Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
% NEED TO FIGURE OUT HOW TO GET v
for i=2:length(k)
    dk(i) = k(i)-k(i-1);
end
dk = dk(2:end);
dk = mean(dk);

xi_sqr = beta.^2.*cumtrapz(waveSpectrum.*abs(Tv_k).^2).*dk;
p_s_ql = exp(-(k_x).^2.*xi_sqr);
end

function [beta,look] = getBeta(metadata,orbit_vec_start,num_orbit_vectors,velocity_value_choice)
    req_atributes = ["slant_range_to_first_pixel","antenna_pointing"];
    meta_beta = filterAttributesNetCDF(metadata.Attributes, req_atributes);
    slant_range = meta_beta(1).Value;
    look = meta_beta(2).Value;
    velocity = getSatVelocity(metadata,num_orbit_vectors,orbit_vec_start,velocity_value_choice);
    beta = slant_range./velocity;
end

function velocity = getSatVelocity(metadata,n,orbit_vec_start,velocity_value_choice)
% Velocity in m/s
updated_attributes = cell(3, n);

% Loop through values from 1 to 17
for i = 1:n    % Construct the attribute strings with the current value for x_vel, y_vel, and z_vel
    x_vel_str = sprintf('Orbit_State_Vectors:orbit_vector%d:x_vel', i);
    y_vel_str = sprintf('Orbit_State_Vectors:orbit_vector%d:y_vel', i);
    z_vel_str = sprintf('Orbit_State_Vectors:orbit_vector%d:z_vel', i);
    
    % Add the attribute strings to the cell array
    updated_attributes{1, i} = x_vel_str;
    updated_attributes{2, i} = y_vel_str;
    updated_attributes{3, i} = z_vel_str;
end

% Convert the cell array to comma-separated strings
x_vel_attributes = strjoin(updated_attributes(1, :), ',');
y_vel_attributes = strjoin(updated_attributes(2, :), ',');
z_vel_attributes = strjoin(updated_attributes(3, :), ',');

% Join the x_vel, y_vel, and z_vel attributes into a single string
all_attributes = strcat(x_vel_attributes, ',', y_vel_attributes, ',', z_vel_attributes);

meta_atr = metadata.Attributes;
%orbit_vec_start = 245;
meta_orb = meta_atr(orbit_vec_start:orbit_vec_start+n*7-1);

start_x_vel = 5;
start_y_vel = 6;
start_z_vel = 7;

x_vel_indices = linspace(start_x_vel,start_x_vel+(n-1)*7,n);
y_vel_indices = linspace(start_y_vel,start_y_vel+(n-1)*7,n);
z_vel_indices = linspace(start_z_vel,start_z_vel+(n-1)*7,n);

meta_orb_x = [meta_orb(x_vel_indices).Value];
meta_orb_y = [meta_orb(y_vel_indices).Value];
meta_orb_z = [meta_orb(z_vel_indices).Value];

velocity_calc = sqrt(meta_orb_x.^2 + meta_orb_y.^2 + meta_orb_z.^2);
%disp(velocity_calc);
switch velocity_value_choice
    case 'max'
        velocity = mean(velocity_calc);
    case 'mean'
        velocity = mean(velocity_calc);
    case 'min'
        velocity = min(velocity_calc);
    otherwise
        error("Incorrect velocity value choice. Make sure that the value is either 'min', 'mean', or 'max'.");
end

end