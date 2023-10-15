function info = readEnviHdr(filename)
% readEnviHdr takes the filepath of the desired HDR file and outputs a
% structure containing all the required information for processing using
% the multibandread function. The values are converted to integers and the
% correct strings are added to replace the values used for 'data types' and
% 'byte order'. The associated structure field names match that of the ENVI
% file as well as an additional fields, called 'size' and 'file'
% size: places the samples, lines, and bands in array of the form: [lines,samples,bands]
% file: stores the file path of the associated .img file

%% Check header file exists
fid = fopen(filename);
if fid<0
    error('%s does not exist. \n', filename); 
end

%% Load whole header file into a string
str_info = fread(fid,'uint8=>char')';
fclose(fid);

%% String manipulation
% Split lines of string into column vector
lines = strsplit(str_info, '\n').';
lines = lines(2:end-1);

%% Fill struc with appropriate values from column vector
info = struct();
for lineIndex=1:length(lines)
    [param, value] = LineFormat(lines{lineIndex});
    info.(param) = value;
end

%% Translate appropriate values to associated strings for multibandread
% Update 'byte-order' field
if isfield(info,'byte_order')
    switch info.byte_order
        case 0
            info.byte_order = 'ieee-le';
        case 1
            info.byte_order = 'ieee-be';
        otherwise
            info.byte_order = 'n';
    end
end

% Update 'data_type' field

if isfield(info,'data_type')
    switch info.data_type
        case 1
            info.data_type = 'int8';
        case 2
            info.data_type = 'int16';
        case 3
            info.data_type = 'int32';
        case 4
            info.data_type = 'float';
        case 12
            info.data_type = 'uint16';
        otherwise
            error(['Unsupported file type number: ', num2str(info.value)])
    end
end

%% Create size array for multibandread
if isfield(info, 'lines') && isfield(info, 'samples') && isfield(info, 'bands')
  info.size = [info.lines, info.samples, info.bands];
end

%% Create file array for multibandread
extIndex = find(filename=='.',1,"first");
if ~isempty(extIndex)
    fileName = strtrim(filename(1:extIndex-1)); % extract file name
    fileName = [fileName '.img']; % add .img extension
    info.file = fileName;
end

end

function [param, value] = LineFormat(line)
% LineFormat takes in a line of the ENVI file and formats the string to
% extract the parameter and value as either a string or integer value,
% where applicable. The parameter names are updated to remove spaces and
% add an underscore character for formatting in a structure
    param = '';
    valueChar = '';
    eqIndex = find(line=='=',1,"first");
    if ~isempty(eqIndex)
        param = strtrim(line(1:eqIndex-1)); % extract parameter name
        param = strrep(param,' ', '_'); % Update field name to remove spaces
        valueChar = strtrim(line(eqIndex+1:end)); % extract value
    end
    valueNum = str2double(valueChar);
    if ~isnan(valueNum) % Check if the value is a number
        value = valueNum;
    else
        value = valueChar;
    end
end
