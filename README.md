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


#### [`report`](./docs/report/)


### [`functions`](./functions/)

#### [`hasselmann`](./functions/hasselmann)

#### [`plotting`](./functions/plotting)

#### [`waveSpectra`](./functions/waveSpectra)


### [`plots`](./plots/)

### [`toolboxes`](./toolboxes/)

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
3. wgrib
  ```sh
  npm install npm@latest -g
  ```

### Installation

## Configuration


## Usage


## License

## Contact

