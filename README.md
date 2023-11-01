# SAR Parameter Extraction Pipeline
This repository contains the MATLAB code that was written for an undergraduate final-year project (EEE4022S) at the University of Cape Town. The aim of the project was to perform an initial investigation into designing a parameter extraction pipeline for SAR images to extract the sea ice characteristics of Antarctica.

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#repo-overview">Repo Overview</a>
      <ul>
        <li><a href="#docs">docs</a></li>
        <li><a href="#functions">functions</a></li>
        <li><a href="#plots">plots</a></li>
        <li><a href="#toolboxes">toolboxes</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

## Repo Overview

A description of the various folders and functions contained within this directory is given below.
### [`docs`](./docs/)


#### [`report`](./docs/JNSRYA006_EEE4022S_RAV2023-03.pdf/)
This is the submitted report for this project and contains an overview of the design and implementation procedure.

#### [`pipeline`](./docs/pipeline.pdf/)
This is a PDF version of the MATLAB live-script. Known as the pipeline.

### [`functions`](./functions/)
This folder contains all functions used in implementing the pipeline. The folder is broken down into three sub-folders, each containing the respective functions for the associated folder name.

#### [`hasselmann`](./functions/hasselmann)
This folder contains all functions used in implementing the [Hasselmann and Hasselmann](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/91JC00302) and its inversion. These functions are implemented in the pipeline and can be individually run.

#### [`plotting`](./functions/plotting)
This folder contains all functions used for generating plots of the respective spectra and validation plots.

#### [`waveSpectra`](./functions/waveSpectra)
This folder contains all functions used to generate the respective spectra used as an input to the Hasselmann and Hasselmann inversion procedure. Data used in generating these spectra can be obtained from [NOAA NCEP](https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/).

### [`toolboxes`](./toolboxes/)
This folder contains all external toolboxes utilised in the development of this pipeline. Namely, MMaps, and a JONSWAP wave generation spectrum, sourced from [Thor Fossen](https://github.com/cybergalactic/MSS).

## Getting Started

### Prerequisites
The following software and MATLAB toolboxes are required to run this pipeline.
1. #### SNAP ESA
SNAP is a tool developed by the European Space Agency (ESA) to process SAR data obtained from Sentinel Satellites. SNAP is used to pre-process data in this pipeline and is essential to install.
To install SNAP please do the following:
- Download [SNAP here](https://step.esa.int/main/download/snap-download/) for your OS distribution
- Choose only the Sentinel Toolboxes installer
- Install SNAP and follow the onscreen instructions. It is only necessary to install the Sentinel-1 Toolbox
2. #### M_Map
M_Map is a mapping package for MATLAB which allows data plots on different world maps.
To use M_Map please do the following:
- Download the zipped package here and
- Add the extracted folder to your MATLAB path
- _Note:_ A copy of the M_Map package is included on this git repo in the following path: [`toolboxes/m_map`](./toolboxes/m_map)
3. #### wgrib2
wgrid2 is a utility package developed by NOAA NCEP and is [available here](https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/).
Installation differs depending on your OS. To install and compile wgrib2 for a different OS distribution use the appropriate links.
- [For Windows](https://ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/Windows10/Installation)
- [For Linux](https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/compile_questions.html)
- [For MacOS](https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/compile_questions.html)

### Installation
In order to install and run this pipeline, the following steps need to be followed.
1. Clone the repo
2. Open MATLAB and navigate to the _sar-parameter-extraction-pipeline_ folder
   
## Configuration
1. Add **all folders and subfolders** to the path

## Usage
1. Open [pipeline.mlx](./functions/pipeline.mlx/)
2. Import pre-processed SAR data. Make sure the instructions in the Live-Script are followed.
3. Run each section. Plots can be toggled on or off by the respective tick boxes.


## License
This repository and all work herein are distributed under an MIT License.

## Contact
For any issues or queries, please contact Ryan via [email](mailto:JNSRYA006@myuct.ac.za).

