import pygrib
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--file_path', type=str)
args = parser.parse_args()

grbs = pygrib.open(args.file_path)
#grbs =  pygrib.open(r'C:\PhD\Radar\ADLP\MRMS\RawData\2020-07-23\PrecipRate_00.00_20200723-234000.grib2')
grbs.seek(0)
grb = grbs[1]

Values = grb.values
# Extract latitudes and longitudes
lats, lons = grb.latlons()
lons=lons-360



# Define latitude and longitude ranges
lat_min, lat_max = 33.97757, 34.22959
lon_min, lon_max = -106.9664, -106.6382


# Find the indices of the latitudes and longitudes within the defined ranges
lat_indices = np.where((lats >= lat_min) & (lats <= lat_max))[0]
lon_indices = np.where((lons >= lon_min) & (lons <= lon_max))[1]

# Use the indices to extract the values within the defined ranges
Radar_data = Values[lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1]
lon = lons[lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1]
lat = lats[lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1]
