function fh = helperFunctions
fh.getCaptureDate = @getCaptureDate;
fh.getNOAAParams = @getNOAAParams;
fh.look = @lookDiscretise;
fh.getLook = @getLook;
fh.getPolarisation = @getPolarisation;
fh.getCaptureTime = @getCaptureTime;
fh.getBeta = @getBeta;
fh.getSatVelocity = @getSatVelocity;
fh.kl = @defineKLook;
fh.omega = @defineOmega;
fh.resize = @resizeToSameSize;
fh.tiltMTF = @tiltMTF;
fh.hydroMTF = @hydrodynamicMTF;
fh.rangeVelocityTF = @rangeVelocityTF;
fh.rarMTF = @rarMTF;
fh.velocityBunchingMTF = @velocityBunchingMTF;
fh.sarImagingMTF = @sarImagingMTF;
end

function [noaaDateStr, noaaHourStr] = getNOAAParams(SARCaptureDate)
    % Extract year, month, day and hour, and minutes
    noaaDate = datetime(SARCaptureDate,"Format","uuuuMMdd");
    noaaHour = datetime(SARCaptureDate,"Format","HH:mm");
    % Convert hour to double
    noaaHour = double(hour(noaaHour));
    % Set hourly intervals of NOAA NCEP Wave data
    allowedHours = [0,6,12,18];
    % Calulate the smallest differene between SAR Capture Time and NOAA
    % NCEP time intervals and get index of closest hour
    hourDiff = abs(allowedHours - noaaHour);
    [~,hourIndex] = min(hourDiff);
    % Select closest hour from allowedHours
    nearestHour = allowedHours(hourIndex);
    % Isolate month and day from capture date
    noaaMonth = num2str(month(noaaDate));
    noaaDay = num2str(day(noaaDate));
    % Check if single digit values (month and day functions on return
    % single character if date or minth value < 10. NOAA NCEP expects every
    % day and month character to be two digits.
    if length(noaaMonth) == 1
        noaaMonth = ['0',noaaMonth];
    end
    if length(noaaDay) == 1
        noaaDay = ['0',noaaDay];
    end      
    % Return strings in correct format for downloading wave data
    noaaDateStr = [num2str(year(noaaDate)),noaaMonth,noaaDay];
    noaaHourStr = num2str(nearestHour);
end

function captureDate = getCaptureDate(metadata)
    reqAttributes = "PROC_TIME";
    metaDate = filterAttributesNetCDF(metadata.Attributes, reqAttributes);
    format = 'dd-MMM-yyyy HH:mm:ss.SSSSSS';
    captureDate = datetime(metaDate(1).Value, 'InputFormat', format);
end

function captureTime = getCaptureTime(metadata)
    reqAttributes = ["first_line_time","last_line_time"];
    meta_time = filterAttributesNetCDF(metadata.Attributes, reqAttributes);
    format = 'dd-MMM-yyyy HH:mm:ss.SSSSSS';
    captureTime(1) = datetime(meta_time(1).Value, 'InputFormat', format);
    captureTime(2) = datetime(meta_time(2).Value, 'InputFormat', format);
    captureTime = mean(captureTime);
end

function polarisation = getPolarisation(metadata)
    reqAttributes = ["mds1_tx_rx_polar","mds2_tx_rx_polar"];
    meta_polar = filterAttributesNetCDF(metadata.Attributes, reqAttributes);
    polarisation = ["",""];
    polarisation(1) = meta_polar(1).Value;
    polarisation(2) = meta_polar(2).Value;    

end

function look = getLook(metadata)
    reqAttributes = "antenna_pointing";
    meta_look = filterAttributesNetCDF(metadata.Attributes, reqAttributes);
    look = meta_look(1).Value;
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

function beta = getBeta(metadata)
    req_atributes = "slant_range_to_first_pixel";
    meta_beta = filterAttributesNetCDF(metadata.Attributes, req_atributes);
    slant_range = meta_beta(1).Value;
    orbit_vector_time = getCaptureTime(metadata);
    velocity = getSatVelocity(metadata,orbit_vector_time);
    beta = slant_range./velocity;
end

function velocity = getSatVelocity(metadata,orbit_vector_time)
% Velocity in m/s
meta_atr = metadata.Attributes;
for i = 1:numel(meta_atr)
    % Check if the 'Name' field matches the field name you want
    if strcmp(meta_atr(i).Name, 'Abstracted_Metadata:Orbit_State_Vectors:orbit_vector1:time')
        % Store the index and break the loop if a match is found
        orbit_vec_start = i;
        break;
    end
end
for i = 1:numel(meta_atr)
    % Check if the 'Name' field matches the field name you want
    if strcmp(meta_atr(i).Name,'Abstracted_Metadata:SRGR_Coefficients:srgr_coef_list_1:zero_doppler_time')
        % Store the index and break the loop if a match is found
        orbit_vec_end = i;
        break;
    end
end
n = (orbit_vec_end - orbit_vec_start)/7;

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

meta_orb = meta_atr(orbit_vec_start:orbit_vec_start+n*7-1);

start_time = 1;
start_x_vel = 5;
start_y_vel = 6;
start_z_vel = 7;

time_indices = linspace(start_time,start_time+(n-1)*7,n);
x_vel_indices = linspace(start_x_vel,start_x_vel+(n-1)*7,n);
y_vel_indices = linspace(start_y_vel,start_y_vel+(n-1)*7,n);
z_vel_indices = linspace(start_z_vel,start_z_vel+(n-1)*7,n);

meta_orb_time = [meta_orb(time_indices).Value];
format = 'dd-MMM-yyyy HH:mm:ss.SSSSSS';
meta_orb_time_dt = datetime(meta_orb_time, 'InputFormat', format);

% Calculate time differences
time_diff = abs(meta_orb_time_dt - orbit_vector_time);
% Find the index of the closest datetime
[~, time_index] = min(time_diff);

meta_orb_x = [meta_orb(x_vel_indices(time_index)).Value];
meta_orb_y = [meta_orb(y_vel_indices(time_index)).Value];
meta_orb_z = [meta_orb(z_vel_indices(time_index)).Value];

velocity = sqrt(meta_orb_x.^2 + meta_orb_y.^2 + meta_orb_z.^2);
%disp(velocity);

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

function TR_k = rarMTF(Tt_k,Th_k)
    TR_k = Tt_k + Th_k;
end

function Tvb_k = velocityBunchingMTF(beta,k_x,Tv_k)
    Tvb_k = -1i.*beta.*k_x.*Tv_k;
end

function TS_k = sarImagingMTF(TR_k,Tvb_k)   
    TS_k    = TR_k + Tvb_k;
end