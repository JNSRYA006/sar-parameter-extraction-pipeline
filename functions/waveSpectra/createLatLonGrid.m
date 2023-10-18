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