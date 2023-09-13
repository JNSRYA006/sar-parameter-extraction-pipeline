

function [all_metadata] = readMetadata(metadata)

% Need to figure out which orbit number is used
% What is frames?
% need to calculate pixel size in m?
% calculate velocity
% heading
% LOOK_SEP??
% Slant range
% angle of incidence?? which value will it be
% everything after n-spec
req_atributes = ["MISSION","orbit_cycle","first_line_time","antenna_pointing","PASS","centre_heading","slant_range_to_first_pixel","centre_lat","centre_lon","total_size"];
filtered_metadata = filterAttributesH5(metadata,req_atributes);
all_metadata = getAttributeValH5(filtered_metadata,req_atributes);

end