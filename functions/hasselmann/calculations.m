%% Calculate Oribital Velocity Covariance
% Set file location of H5 file
fileLoc = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt.h5';
testOut = getMetadataH5(fileLoc, 'abstracted');
testOutAtr = testOut.Attributes;

req_atributes = ["antenna_pointing","incidence_near","incidence_far"];
test_output = filterAttributesH5(testOutAtr,req_atributes);
test_meta_val = getAttributeValH5(test_output,req_atributes);
%%
look_index = strcmp({test_meta_val.Name},'antenna_pointing');
test_meta_val(1).Value;
look = lookDiscretise(test_meta_val(1).Value);
incidence_near_index = strcmp({test_meta_val.Name},'incidence_near');
incidence_near = test_meta_val(incidence_near_index).Value;
incidence_far_index = strcmp({test_meta_val.Name},'incidence_far');
incidence_far = test_meta_val(incidence_far_index).Value;
th = extrapolateIncidence(incidence_near,incidence_far,100);

F_k = 0;
fv_k = orbitalVelocityCovariance(F_k,th,look,heading);

function lookVal = lookDiscretise(look)
% Returns 0 for right look, 1 for left look

switch look

    case 'right'
        lookVal = 0;
    case 'left'
        lookVal = 1;
    otherwise
        error('Invalid input: "%s"', look);
end
end

function incidence = extrapolateIncidence(incidence_near, incidence_far, num_pixels)
% Creates a linscape array of the incidence angle as it changes through the
% vertical number of pixels
if (ischar(incidence_near))
    incidence_near = str2double(incidence_near);
end

if (ischar(incidence_far))
    incidence_far = str2double(incidence_far);
end

incidence = linspace(incidence_near,incidence_far,num_pixels);

end