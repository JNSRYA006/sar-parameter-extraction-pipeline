function [p_s_ql,beta,xi_sqr] = quasilinearCoeff(k,k_y,k_x,waveSpectrum,SARmetadata,th)

func = helperFunctions;

beta = func.getBeta(SARmetadata);
look = func.getLook(SARmetadata);

look = func.look(look);
%th = ncread(SARmetadata.Filename,"Incidence_Angle");
%th = func.incidence(inc_near,inc_far,num_pixels);
%th = extrapolateIncidence(inc_near,inc_far,num_pixels); % Need to calculate based on velocity and look time
% Need to check resizing is correct
%k_new = func.resize(k,th(1,:));
k_new = k;
%k_y = func.resize(k_y,th(1,:));
%k_x = func.resize(k_x,th(1,:));
k_l = func.kl(look,k_y);
omega = func.omega(k_new);
%% Could be wrong
% Create a logical mask for rows with NaN entries
nanRows = any(isnan(waveSpectrum), 2);

% Remove rows with NaN entries
waveSpectrum = waveSpectrum(2:end,2:end);

waveSpectrum = func.resize(waveSpectrum,th);
Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tv_k = func.resize(Tv_k(:,2:end),waveSpectrum);
% NEED TO FIGURE OUT HOW TO GET v
for i=2:length(k)
    dk(i) = k(i)-k(i-1);
end
dk = dk(2:end);
dk = mean(dk);

for i=2:length(k_x)
    dk_x(i) = k_x(i)-k_x(i-1);
end
dk_x = dk_x(2:end);
dk_x = abs(mean(dk_x));

for i=2:length(k_y)
    dk_y(i) = k_y(i)-k_y(i-1);
end
dk_y = dk_y(2:end);
dk_y = abs(mean(dk_y));

xi_sqr = beta.^2.*cumtrapz(cumtrapz(waveSpectrum.*abs(Tv_k).^2).*dk_x).*dk_y;
xi_sqr = func.resize(xi_sqr(:,2:end),waveSpectrum);
xi_sqr = beta.^2.*trapz(trapz(waveSpectrum.*abs(Tv_k).^2).*dk_x).*dk_y;
xi = sqrt(xi_sqr);
test_k_neg = -(k_x).^2;
k_x_cutoff = xi^(-1);

% Cutoff k_x
% Initialize the index variable
index = 0;
% Iterate through the array to find the first index greater than the threshold
for i = 1:length(k_x)
    if k_x(i) > k_x_cutoff  
        index = i;
        break;  % Exit the loop when the first index is found
    end
end

p_s_ql = exp(-(k_x).^2*xi_sqr);
if ~index == 0
    p_s_ql = func.resize(p_s_ql(1:index),waveSpectrum(1,:));
end
end

function beta = getBeta(metadata,orbit_vec_start,num_orbit_vectors,velocity_value_choice)
    req_atributes = ["slant_range_to_first_pixel"];
    meta_beta = filterAttributesNetCDF(metadata.Attributes, req_atributes);
    slant_range = meta_beta(1).Value;
    %look = meta_beta(2).Value;
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