# SAR Parameter Extraction Pipeline
This directory contains the MATLAB code that was written for an undergraduate final-year project (EEE4022S) at the University of Cape Town. The aim of the project was to perform an initial investigation into designing a parameter extraction pipeline for SAR images to extract the sea ice characteristics of Antarctica.

This folder contains all functions used in implementing the Hasselmann procedure. The respective co- and auto-covariance functions, as well as the spectral expansions, are shown in this directory.

##### [`helperFunctions`](./helperFunctions.m)
This function structure contains all utility functions used to calculate the required Modulation Transfer Functions (MTF)s in the implementation. All the additional metadata functions are contained herein. The functions are called using the following command:

```
func = helperFunctions;
look = func.getLook(metadata) % Example call of a function contained in the helperFunctions function structure
```
##### [`generateSARSpectrumOceanWaves`](./generateSARSpectrumOceanWaves.m)
This function generates the equivalent SAR spectrum of a given ocean wave spectrum according to the Hasselmann procedure.

##### [`inversion`](./inversion.m)
This function acts as the inversion of the Hasselmann procedure and iterates the cost function for a user-defined number of iterations.
