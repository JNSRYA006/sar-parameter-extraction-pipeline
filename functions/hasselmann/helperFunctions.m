function fh = helperFunctions
fh.incidence = @extrapolateIncidence;
fh.look = @lookDiscretise;
fh.kl = @defineKLook;
fh.omega = @defineOmega;
fh.resize = @resizeToSameSize;
fh.tiltMTF = @tiltMTF;
fh.hydroMTF = @hydrodynamicMTF;
fh.rangeVelocityTF = @rangeVelocityTF;
end

function incidence = extrapolateIncidence(incidence_near,incidence_far,num_pixels)
% Creates a linscape array of the incidence angle as it changes through the
% vertical number of pixels in degrees

    if (ischar(incidence_near))
        incidence_near = str2double(incidence_near);
    end

    if (ischar(incidence_far))
        incidence_far = str2double(incidence_far);
    end

    incidence = linspace(incidence_near,incidence_far,num_pixels);

end

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

function k_l = defineKLook (lookVal,k_y)
% Sets k_l depending on whether the satellite is left or right looking
if ~(lookVal)
    k_l = -k_y;
    return;
end
k_l = k_y;
end

function omega = defineOmega(k)
    g = 9.81;
    omega = abs(sqrt(g.*k));
end

function resizedMatrix = resizeToSameSize(matrixToResize, referenceMatrix)
    % Get the size of the reference matrix
    [rows, cols] = size(referenceMatrix);
    
    % Resize the matrix to match the size of the reference matrix
    resizedMatrix = imresize(matrixToResize, [rows, cols]);
end

function Tt_k = tiltMTF(polarisation,k_l,incidenceAngle)
status = 0;
for i = 1:length(polarisation)
    switch polarisation(i)
        case "HH"
            Tt_k = 4.*1i.*k_l.*cotd(incidenceAngle).*(1+sind(incidenceAngle).^2).^(-1);
            status = 0;
        case "VV"
            Tt_k = 8.*1i.*k_l.*(sind(2.*incidenceAngle)).^(-1);
            status = 0;
        otherwise
            status = 1;
    end
end
if(status)
    error("SAR data are in neither VV or HH polarisations and Tilt MTF cannot be computed");
end
end

function Th_k = hydrodynamicMTF(omega,mu,k,k_y)
Th_k = (omega-1i.*mu)./(omega.^2+mu.^2).*(4.5).*k.*omega.*((k_y.^2)./(k.^2));
end


function Tv_k = rangeVelocityTF(omega,incidenceAngle,k_l,k)

Tv_k = -omega.*(sind(incidenceAngle).*k_l./abs(k) + 1i.*cosd(incidenceAngle));

end