function metadata = filterAttributesNetCDF(metadata_struc, attribute)
% filterAttributesNetCDF takes in a structure containing metadata (obtained
% using meta_nc = ncinfo(filepath,'metadata');), as well as a list of strings of the
% desired attributes to filter. A list of available attribute types are
% available in the Appendix of the report. 
% The function outputs a 1xlength(attributes) structure with all the
% original metadata columns preserved
prefix = 'Abstracted_Metadata:';
attributes = split(attribute, ',');
num_attributes = length(attributes);
metadata = struct([]);

%% Update attribute array to include prefix for abstracted metadata
attribute_update = strings(num_attributes);
for i=1:num_attributes
    attribute_update(i) = prefix + attribute(i);
end

%% Extract input attributes and filter out columns
for i = 1:num_attributes
    index = strcmp({metadata_struc.Name}, attribute_update(i));
    metadata = [metadata, metadata_struc(index)];
end

%% Remove empty rows and prefix
%metadata.Name = attributes;
for i = 1:num_attributes
    % Remove the common prefix by indexing the string
    metadata(i).Name = attributes(i);
end

metadata = removeEmptyRows(metadata);
end

function cleaned_structure = removeEmptyRows(structure)
%% Remove any empty rows in a structure
isEmptyRow = false(size(structure));

for j = 1:numel(structure)
    isEmptyRow(j) = all(structfun(@isempty,structure(j)));
end
cleaned_structure = structure(~isEmptyRow);
end