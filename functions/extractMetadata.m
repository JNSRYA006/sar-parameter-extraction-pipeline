fileLoc = 'D:\UCT\EEE4022S\Data\CPT\subset_2_cpt.h5';
testOut = getMetadataH5(fileLoc, 'abstracted');
testOutAtr = testOut.Attributes;

%% Testing
attribute_names = ["MISSION","SWATH", "BEAMS","PRODUCT", "orbit_cycle", "ABS_ORBIT"];
attribute_names2 = ["MISSION","SWATH", "BEAMS", "ABS_ORBIT"];
test_output = filterAttributesH5(testOutAtr,attribute_names);

test_meta_val = getAttributeValH5(test_output,attribute_names2);

% pro F_read_parameters, SAR_log_file = SAR_log_file, $
%                         WAVE_log_file = sea_log_file, $
%                         SENSOR, $
%                         ORBIT, $ ; ORBIT NUMBER
%                         FRAMES, $ ; FRAMES NUMBER
%                         DAY, $ ; SAR IMAGE ACQUISITION DAY
%                         MONTH, $ ; SAR IMAGE ACQUISITION MONTH
%                         YEAR, $ ; SAR IMAGE ACQUISITION YEAR
%                         LOOK_SAR_, $ ; 'LEFT' or 'RIGHT'
%                         PASS, $ ; 'ASC' or 'DESC'
%                         PIX_SIZE_, $
%                         VEL_SAR, $
%                         HEADING_DEG, $
%                         LOOK_SEP, $
%                         SLANT_RANGE, $
%                         TETA_DEG, $ ; angle of incidence in deg
%                         LAT0_DEG, $ ; tile center latitude in deg
%                         LON0_DEG, $ ; longitude inside tiles in deg
%                         N_SPEC, $ ; number of pixels/size of the SAR spectrum
%                         NCOL_INI, $
%                         NRIG_INI, $
%                         WIN_NUMBER, $
%                         DIR_INP_SAR, $ ; directory where the CO and CROSS SAR spectrum reside
%                         FILE_CROSS_SPEC_SAR, $
%                         FILE_CO_SPEC_SAR, $
%                         DIR_INP_WAVE, $ ; directory where the wave spectra reside
%                         FILE_FULL_WAVE_SPEC