# SAR Parameter Extraction Pipeline
This directory contains the MATLAB code that was written for an undergraduate final-year project (EEE4022S) at the University of Cape Town. The aim of the project was to perform an initial investigation into designing a parameter extraction pipeline for SAR images to extract the sea ice characteristics of Antarctica.

This folder contains all functions used in the wave spectrum generation step of the pipeline. These wave spectra serve as the first-guess input to the Hasselmann procedure.

##### [`createLatLonGrid`](./createLatLonGrid.m)
This function creates a 1x3 grid of latitude and longitude values based on the resolution of wave data available, as well as the location of the point of interest. The grid values are the closest 3 latitude and longitude values to the point of interest.

##### [`generate2DWaveSpectrum`](./generate2DWaveSpectrum.m)
This function generates a two-dimensional frequency-direction wave spectrum based on the one-dimensional JONSWAP spectrum, and the direction distribution function generated using the [`generateDirectionalDistribution`](./generateDirectionalDistribution.m) function.

##### [`generateDirectionalDistribution`](./generateDirectionalDistribution.m)
This function creates a directional distribution function based on the input wave spectrum.

##### `generateJONSWAP`
The [`generateMultipleJONSWAP`](./generateMultipleJONSWAP.m) and [`generateSingleJONSWAP`](./generateSingleJONSWAP.m) functions generate a JONSWAP model of the input wave spectra for multiple and single latitude and longitude points respectively.

##### [`getSubsetWaveVals`](./getSubsetWaveVals.m)
This function takes a subset of the downloaded wave values based on the location of interest. Start and end latitude and longitude coordinates are required, and these can be generated using the [`createLatLonGrid`](./createLatLonGrid.m) function.

##### [`waveNumberSpectrum`](./waveNumberSpectrum.m)
This function creates the associated wave-number spectrum from the frequency-direction spectrum. This wave-number spectrum is used as the input to the Hasselmann generation of an equivalent SAR spectrum of ocean waves.

