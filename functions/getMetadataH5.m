function metadata = getMetadataH5(filename, type)
% getMetadataH5 takes in the .h5 file and the type of metadata desired and
% returns this metadata in a structure. The allowable types are:
% 'abstracted': Abstracted_Metadata
% 'original': Original_Product_Metadata
% 'processing': Processing_Graph
% 'history': history
% 'all': all of the above types
% In order to get the attributes, please use the getAttributesH5 function

%% Check file exists
fid = fopen(filename);
if fid<0
    error('%s does not exist. \n', filename); 
end

%% Read in file
hinfo = h5info(filename);
hinfo_meta = hinfo.Groups(3);
switch type
    case 'abstracted'
            metadata = hinfo_meta.Groups(1);
    case 'original'
            metadata = hinfo_meta.Groups(2);
    case 'processing'
            metadata = hinfo_meta.Groups(3);
    case 'history'
            metadata = hinfo_meta.Groups(4);
    case 'all'
            metadata = hinfo_meta;
    otherwise
            metadata = hinfo_meta;
            error(['%s is not a valid metadata type. \n ' ...
                'Returning all metadata types as a structure... \n'], type); 
end
end