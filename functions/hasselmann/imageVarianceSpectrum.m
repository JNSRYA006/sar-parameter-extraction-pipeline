function p_s_k = imageVarianceSpectrum(k,k_x,k_y,waveSpectrum,waveSpectrum_inv,SARmetadata,th)
% Equation 26 in Hasselmann
%% Get required metadata
func = helperFunctions;
look = func.getLook(SARmetadata);
% Update metadata function to general get look
%th = ncread(SARmetadata.Filename,"Incidence_Angle");

look = func.look(look);
polarisation = func.getPolarisation(SARmetadata);
beta = getBeta(SARmetadata,245,16,'max');

% Need to check resizing is correct
k_new = func.resize(k,th(1,:));
k_y = func.resize(k_y,th(1,:));
k_x = func.resize(k_x,th(1,:));
k_l = func.kl(look,k_y);
omega = func.omega(k_new);
mu=0.5;

Tt_k = func.tiltMTF(polarisation,k_l,th);
Th_k = func.hydroMTF(omega,mu,k_new,k_y);
Tt_k_inv = func.tiltMTF(polarisation,k_l,th);
Th_k_inv = func.hydroMTF(omega,mu,-k_new,k_y);
Tr_k = func.rarMTF(Tt_k,Th_k);
Tr_k_inv = func.rarMTF(Tt_k_inv,Th_k_inv);

Tv_k = func.rangeVelocityTF(omega,th,k_l,k_new);
Tv_k_inv = func.rangeVelocityTF(omega,th,k_l,-k_new);
Tvb_k = func.velocityBunchingMTF(beta,k_x,Tv_k);
Tvb_k_inv = func.velocityBunchingMTF(beta,k_x,Tv_k_inv);

Ts_k = func.sarImagingMTF(Tr_k,Tvb_k);
Ts_k_inv = func.sarImagingMTF(Tr_k_inv,Tvb_k_inv);

% Remove NaN values manually
waveSpectrum = waveSpectrum(2:end,2:end);

% Resize the matrix to 2001x2001
waveSpectrum = func.resize(waveSpectrum,th);

% Remove NaN values manually
waveSpectrum_inv = waveSpectrum_inv(2:end,2:end);

% Resize the matrix to 2001x2001
waveSpectrum_inv = func.resize(waveSpectrum_inv,th);
    
p_s_k = abs(Ts_k).^2.*(waveSpectrum./2) + abs(Ts_k_inv).^2.*(waveSpectrum_inv./2);
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