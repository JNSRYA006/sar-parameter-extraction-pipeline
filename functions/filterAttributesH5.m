function metadata = filterAttributesH5(metadata_struc, attribute)
% filterAttributesH5 takes in a structure containing metadata (obtained
% using the getMetadataH5 function), as well as a list of strings of the
% desired attributes to filter. A list of available attribute types are
% available here: 
% The function outputs a 1xlength(attributes) structure with all the
% original metadata columns preserved
attributes = split(attribute, ',');
num_attributes = length(attributes);
metadata = struct([]);

%% Check if attribute is part of the structure
% if ~isfield(metadata_struc, attribute)
%         error('Attribute, %s, does not exist in the structure, %s. \n', attribute,metadata_struc); 
% end
%test_val = metadata_struc(metadata_struc.Name == 'MISSION');

%% Extract input attributes and filter out columns
for i = 1:num_attributes
    index = strcmp({metadata_struc.Name}, attributes(i));
    metadata = [metadata, metadata_struc(index)];
end

%% Remove empty rows
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