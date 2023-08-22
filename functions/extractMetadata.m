fileLoc = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt.h5';
testOut = getMetadataH5(fileLoc, 'abstracted');
testOutAtr = testOut.Attributes;

%% Testing
attribute_names = ["MISSION","SWATH", "BEAMS","PRODUCT", "orbit_cycle", "ABS_ORBIT"];
attribute_names2 = ["MISSION","SWATH", "BEAMS", "ABS_ORBIT"];
test_output = filterAttributesH5(testOutAtr,attribute_names);

test_meta_val = getAttributeValH5(test_output,attribute_names2);
