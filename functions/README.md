# SAR Parameter Extraction Pipeline
This directory contains the MATLAB code that was written for an undergraduate final-year project (EEE4022S) at the University of Cape Town. The aim of the project was to perform an initial investigation into designing a parameter extraction pipeline for SAR images to extract the sea ice characteristics of Antarctica.

This folder contains all functions used in implementing the pipeline. The folder is broken down into three sub-folders, each containing the respective functions for the associated folder name. A description of the various folders and functions contained within this directory is given below.

#### [`/hasselmann`](./hasselmann)
This folder contains all functions used in implementing the [Hasselmann and Hasselmann](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/91JC00302) procedure and its inversion. These functions are implemented in the pipeline and can be individually run.

#### [`/plotting`](./plotting)
This folder contains all functions used for generating plots of the respective spectra and validation plots.

#### [`/waveSpectra`](./waveSpectra)
This folder contains all functions used to generate the respective spectra used as an input to the Hasselmann and Hasselmann inversion procedure. Data used in generating these spectra can be obtained from [NOAA NCEP](https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/).

##### [`downloadNOAAWaveFile`](./downloadNOAAWaveFile.m)
This function downloads the associated [NOAA NCEP](https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/) dataset that is related to the provided SAR dataset.

##### [`getGribStruct`](./getGribStruct.m)
This function converts the downloaded .grib2 file from [NOAA NCEP](https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/) into the desired .NetCDF4 format.

##### [`pipeline`](./pipeline.mlx)
This function acts as the pipeline for this project. It contains all required function calls to generate best-fit wave parameters from a first-guess wave spectrum and SAR spectrum.
