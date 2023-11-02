# SAR Parameter Extraction Pipeline
This directory contains MATLAB toolboxes, and external functions used for an undergraduate final-year project (EEE4022S) at the University of Cape Town. The aim of the project was to perform an initial investigation into designing a parameter extraction pipeline for SAR images to extract the sea ice characteristics of Antarctica.

This folder contains all toolboxes required in the implementation of the pipeline.

#### [`MMap`](https://github.com/JNSRYA006/sar-parameter-extraction-pipeline/tree/main/toolboxes/m_map)
M_Map is a mapping package for MATLAB which allows data plots on different world maps.
To use this toolbox, this folder needs to be added to the MATLAB path.

#### [`MATLAB2Tikz`](https://github.com/JNSRYA006/sar-parameter-extraction-pipeline/tree/main/toolboxes/matlab2tikz-master/matlab2tikz-master)
MATLAB2Tikz is a package which allows plots generated in MATLAB to be exported as a .tex file for high-quality plotting in LaTeX.

##### [`wavespec`](./wavespec.m)
This function was created by [Thor Fossen](https://github.com/cybergalactic/MSS) and allows multiple wave spectrum models to be generated. The model of interest was the JONSWAP wave model.
