%% Testing find closest lat and long
latitude = 35.15;     % Input latitude
longitude = -100.7; % Input longitude
resolution = 0.25;   % Resolution in degrees

[closest_lat, closest_lon] = findClosestLatLon(latitude, longitude, resolution);

fprintf('Closest Latitude: %.2f\n', closest_lat);
fprintf('Closest Longitude: %.2f\n', closest_lon);

grid_size = 1;       % Size of the grid around the desired coordinates

[grid_lat, grid_lon] = createLatLonGrid(latitude, longitude, resolution);

%% Calculate orbital velocity
% Time average over period during which scattering element is viewed by SAR
% Initialize an empty cell array to store the updated attributes
% Initialize an empty cell array to store the updated attributes
filepath = "D:\UCT\EEE4022S\Data\CPT\small_subset.nc";
nc_file = ncinfo(filepath);
metadata = ncinfo(filepath,'metadata');
func = helperFunctions;
%% Getting incidence angle per pixel
% Only northern hemisphere has been exported manually for pixel-by pixel
% values
incidence_1 = ncread(filepath,'Incidence_Angle');
incidence_2 = ncread(filepath,'incident_angle');
incidence_2_linspace = double(func.resize(incidence_2,incidence_1));
% Plot values
figure;
imagesc(incidence_1);
title('Exported Incidence Angle Tie-Point Grid from SNAP ESA');
ylabel('Latitude (pixels)');
xlabel('Longitude (pixels)')
cb = colorbar;  % create and label the colorbar
cb.Label.String = 'Incidence Angle, \theta (degrees)';

figure;
imagesc(incidence_2);
title('Incidence Angle Tie-Point Grid from SNAP ESA');
ylabel('Latitude (pixels)');
xlabel('Longitude (pixels)')
cb = colorbar;  % create and label the colorbar
cb.Label.String = 'Incidence Angle, \theta (degrees)';

figure;
imagesc(incidence_2_linspace);
title('Linspace of Incidence Angle Tie-Point Grid from SNAP ESA');
ylabel('Latitude (pixels)');
xlabel('Longitude (pixels)')
cb = colorbar;  % create and label the colorbar
cb.Label.String = 'Incidence Angle, \theta (degrees)';
%%
figure;
subplot(1,3,1);
imagesc(incidence_1);
title('Exported Incidence Angle Tie-Point Grid from SNAP ESA');
ylabel('Latitude (pixels)');
xlabel('Longitude (pixels)')
cb = colorbar;  % create and label the colorbar
cb.Label.String = 'Incidence Angle, \theta (degrees)';
subplot(1,3,2);
imagesc(incidence_2_linspace);
title('Linspace of Incidence Angle Tie-Point Grid from SNAP ESA');
ylabel('Latitude (pixels)');
xlabel('Longitude (pixels)')
cb = colorbar;  % create and label the colorbar
cb.Label.String = 'Incidence Angle, \theta (degrees)';
subplot(1,3,3)
difference = abs(incidence_1 - incidence_2_linspace);
imagesc(difference);
title('Difference between Exported Tie-Point Grid and Linspace of Incidence Angle');
ylabel('Latitude (pixels)');
xlabel('Longitude (pixels)')
cb = colorbar;  % create and label the colorbar
cb.Label.String = 'Difference in Incidence Angle, \theta (degrees)';
%% Getting orbit velocities
updated_attributes = cell(3, 17);

% Loop through values from 1 to 17
for i = 1:17
    % Construct the attribute strings with the current value for x_vel, y_vel, and z_vel
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

% Display the resulting string
disp(all_attributes);

% Initialize a 51x1 string array
stringArray = strings(51, 1);

string_attributes = string(all_attributes);
%%
%req_atributes = ['Orbit_State_Vectors:orbit_vector1:x_vel','Orbit_State_Vectors:orbit_vector1:y_vel',"Orbit_State_Vectors:orbit_vector1:z_vel"];
%meta_orb = filterAttributesNetCDF(metadata.Attributes, string_attributes);
meta_atr = metadata.Attributes;
meta_orb = meta_atr(245:363);

x_vel_indices = linspace(5,110,16);
y_vel_indices = linspace(6,111,16);
z_vel_indices = linspace(7,112,16);

x_pos_indices = linspace(2,107,16);
y_pos_indices = linspace(3,108,16);
z_pos_indices = linspace(4,109,16);

t_indices = linspace(1,106,16);

meta_orb_vel_x = [meta_orb(x_vel_indices).Value];
meta_orb_vel_y = [meta_orb(y_vel_indices).Value];
meta_orb_vel_z = [meta_orb(z_vel_indices).Value];

meta_orb_pos_x = [meta_orb(x_pos_indices).Value];
meta_orb_pos_y = [meta_orb(y_pos_indices).Value];
meta_orb_pos_z = [meta_orb(z_pos_indices).Value];

% Define start and end datetime values
% Brazil (descending.nc)
% startDateTime = datetime('2023-09-27 09:19:38', 'Format', 'yyyy-MM-dd HH:mm:ss');
% endDateTime = datetime('2023-09-27 09:22:08', 'Format', 'yyyy-MM-dd HH:mm:ss');
% Small_subset.nc
startDateTime = datetime('2023-07-28 17:33:19', 'Format', 'yyyy-MM-dd HH:mm:ss');
endDateTime = datetime('2023-07-28 17:35:59', 'Format', 'yyyy-MM-dd HH:mm:ss');

% Define the number of points you want in the linspace
numPoints = 16;

% Create a linspace of datetime values
meta_orb_time = linspace(startDateTime, endDateTime, numPoints);

velocity_calc = sqrt(meta_orb_vel_x.^2 + meta_orb_vel_y.^2 + meta_orb_vel_z.^2);
disp(velocity_calc);
avg_velocity = mean(velocity_calc); %m/s

%% Plots for orbits
%% Position
% x_position
figure;
plot(meta_orb_time,meta_orb_pos_x);
grid on
xlabel('Time')
ylabel('x position (m)')

%y_position
figure;
plot(meta_orb_time,meta_orb_pos_y);
grid on
xlabel('Time')
ylabel('y position (m)')

%z_position
figure;
plot(meta_orb_time,meta_orb_pos_z);
grid on
xlabel('Time')
ylabel('z position (m)')

%% Velocity
figure;
subplot(1,3,1);
plot(meta_orb_time,meta_orb_vel_x);
xlabel('Time')
ylabel('$x$ velocity (m/s)','interpreter','latex')
subplot(1,3,2)
plot(meta_orb_time,meta_orb_vel_y);
xlabel('Time')
ylabel('$y$ velocity (m/s)','interpreter','latex')
%title('Velcoity components at 6.8266$^\circ$N,52.3837$^\circ$W','interpreter','latex') % For over brazil (descending.nc)
title('Velcoity components at 34.8951$^\circ$S,16.6894$^\circ$E','interpreter','latex') % For over CPT (small_subset.nc)
%title('Velcoity components at 37.80621$^\circ$N,10.7693$^\circ$W','interpreter','latex') % For over Northern hem (north_hemisphere.nc)
subplot(1,3,3)
plot(meta_orb_time,meta_orb_vel_z);
xlabel('Time')
ylabel('$z$ velocity (m/s)','interpreter','latex')

%% Fucntions
function [closest_lat, closest_lon] = findClosestLatLon(latitude, longitude, resolution)
    % Define the latitude and longitude grid based on the given resolution
    lat_grid = -90:resolution:90;
    lon_grid = -180:resolution:180;

    % Find the closest latitude value
    [~, lat_index] = min(abs(lat_grid - latitude));
    closest_lat = lat_grid(lat_index);

    % Find the closest longitude value
    [~, lon_index] = min(abs(lon_grid - longitude));
    closest_lon = lon_grid(lon_index);
end

function [lat_range, lon_range] = createLatLonGrid(latitude, longitude, resolution)
    % Define the latitude and longitude grid based on the given resolution
    lat_grid = -90:resolution:90;
    lon_grid = -180:resolution:180;
    grid_size = 1;

    % Find the closest latitude and longitude values
    [~, lat_index] = min(abs(lat_grid - latitude));
    [~, lon_index] = min(abs(lon_grid - longitude));
    
    % Create arrays of latitude and longitude coordinates around the desired coordinates
    lat_range = lat_grid(lat_index-grid_size:lat_index+grid_size);
    lon_range = lon_grid(lon_index-grid_size:lon_index+grid_size);
end