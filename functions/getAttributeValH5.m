function meta_val = getAttributeValH5(metadata_struc,attribute)
% getAttributeValH5 takes in a structure containing metadata (obtained using 
% the getMetadataH5 function, and filtered using the filterAttributesH5
% function), as well as a list of strings of the desired attributes to 
% filter.
% The function outputs a 1xlength(attributes) structure with just the Name
% and Value
attributes = split(attribute, ',');
num_attributes = length(attributes);
meta_val = struct([]);

for i = 1:num_attributes
    index = strcmp({metadata_struc.Name}, attributes(i));
    meta_val(index).Name = metadata_struc(index).Name;
    meta_val(index).Value = metadata_struc(index).Value;
end

%% Remove empty rows
meta_val= removeEmptyRows(meta_val);
end


function cleaned_structure = removeEmptyRows(structure)
%% Remove any empty rows in a structure
isEmptyRow = false(size(structure));

for j = 1:numel(structure)
    isEmptyRow(j) = all(structfun(@isempty,structure(j)));
end
cleaned_structure = structure(~isEmptyRow);
end