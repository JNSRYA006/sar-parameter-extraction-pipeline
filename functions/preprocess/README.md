# SAR Parameter Extraction Pipeline
This directory contains the MATLAB code that was written for an undergraduate final-year project (EEE4022S) at the University of Cape Town. The aim of the project was to perform an initial investigation into designing a parameter extraction pipeline for SAR images to extract the sea ice characteristics of Antarctica.

This folder contains all functions used in the pre-processing step of the pipeline.

##### [`annotate512Transect`](./annotate512Transect.m)
This function can be used to annotate the location of the full SAR scene where the taken transects are located. Text can be added, as well as a background colour behind the text.

##### [`filterAttributesNetCDF`](./filterAttributesNetCDF.m)
This function generates a 1x1 structure containing all the filtered metadata values from the original metadata structure obtained from Sentinel-1A data. The input to this function is a list of strings, which can be found in the [`report`](https://github.com/JNSRYA006/sar-parameter-extraction-pipeline/blob/main/docs/JNSRYA006_EEE4022S_RAV2023-03.pdf).
