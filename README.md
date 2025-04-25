### GLOB

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

GLOB is a multiface solar irradiance instrument composed of 25 pyranometers (Apogee SP-710) plus one looking at the ground. This project is part of my PhD research at The University Centre in Svalbard and is funded by The Arctic Field Grant (Grant nr xxx).

## Table of Contents

- [Features](#features)
- [Pyranometer Orientations](#pyranometer-orientations)
- [Scripts](#scripts)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)

## Features

- Measures solar irradiance from multiple orientations.
- Includes scripts for data analysis, irradiance calculations, and visualization.
- Utilizes the Faiman et al. (1992) method for beam and diffuse irradiance estimation.
- Works with NetCDF data, which can be requested by contacting Arthur G.

## Pyranometer Orientations

The following table describes the variable names in the data files together with the different orientations of the 26 pyranometers:

| Name     | Var_name     | Azimuth | Tilt | Comments                |
|----------|--------------|---------|------|-------------------------|
| GHI      | SIR_Avg[1]   | 0       | 0    | horizontal plane GHI    |
| SIR135S  | SIR_Avg[2]   | 180     | 45   |                         |
| SIR90S   | SIR_Avg[3]   | 180     | 90   |                         |
| SIR45S   | SIR_Avg[4]   | 180     | 135  |                         |
| SIR135SW | SIR_Avg[5]   | 225     | 45   |                         |
| SIR90SW  | SIR_Avg[6]   | 225     | 90   |                         |
| SIR45SW  | SIR_Avg[7]   | 225     | 135  |                         |
| SIR135W  | SIR_Avg[8]   | 270     | 45   |                         |
| SIR90W   | SIR_Avg[9]   | 270     | 90   |                         |
| SIR45W   | SIR_Avg[10]  | 270     | 135  |                         |
| SIR135NW | SIR_Avg[11]  | 315     | 45   |                         |
| SIR90NW  | SIR_Avg[12]  | 315     | 90   |                         |
| SIR45NW  | SIR_Avg[13]  | 315     | 135  |                         |
| SIR135N  | SIR_Avg[14]  | 0       | 45   |                         |
| SIR90N   | SIR_Avg[15]  | 0       | 90   |                         |
| SIR45N   | SIR_Avg[16]  | 0       | 135  |                         |
| SIR135NE | SIR_Avg[17]  | 45      | 45   |                         |
| SIR90NE  | SIR_Avg[18]  | 45      | 90   |                         |
| SIR45NE  | SIR_Avg[19]  | 45      | 135  |                         |
| SIR135E  | SIR_Avg[20]  | 90      | 45   |                         |
| SIR90E   | SIR_Avg[21]  | 90      | 90   |                         |
| SIR45E   | SIR_Avg[22]  | 90      | 135  |                         |
| SIR135SE | SIR_Avg[23]  | 135     | 45   |                         |
| SIR90SE  | SIR_Avg[24]  | 135     | 90   |                         |
| SIR45SE  | SIR_Avg[25]  | 135     | 135  |                         |
| REF      | Kglob        | 0       | 180  | reflected GHI           |


## Scripts

Sorted by order of execution.

### `PROCESS_GLOB_csvdata.py`
Processes raw GLOB data from CSV files (data not provided), standardizes naming conventions, converts timestamps to UTC, and writes processed data to new CSV files called "GLOB_data_30sec_*.dat".

### `CREATE_ncdf_GLOB_data_LYR.py` and `CREATE_ncdf_GLOB_data_NYA.py`
Processes GLOB data from CSV files (see "\Data_GLOB\GLOB_data_30sec_*.dat"), calculates solar angles, computes surface albedo, and converts data to NetCDF format with metadata called "GLOB_data_5min_*.nc".

### `GLOB_estimations_Faiman_LYR.py` and `GLOB_estimations_Faiman_NYA.py`
Calculates beam and diffuse irradiance for GLOB data in Longyearbyen/Ny-Ålesund using Faiman et al. (1992) method. Processes NetCDF data (see "\Data_GLOB\GLOB_data_5min_*.nc") and save the output in CSV files, generates plots, and compares results with reference data.

### `CALC_irradiance_multiorient.py`
Calculates solar irradiance for multiple orientations using Faiman et al. (1992) method. Loads data from CSV files generated with `GLOB_estimations_Faiman_*.py` in a folder called "\MultiOrientations_Irradiance_Data_LYR", computes angles, filters negative values, and saves results to CSV.

### `PLOT_irradiance_multiorient_polar_plot.py`
Generates visualizations, including polar heatmaps and time series plots, to analyze irradiance components from GLOB data.

### `glob_functions_Faiman.py`
Provides functions for calculating solar angles, irradiance components, and geometry coefficients using the Faiman et al. (1992) method. Includes parallel processing for efficient computation.

## Installation

To set up the project locally, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/GLOB.git
    cd GLOB
    ```

2. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

0/ `PROCESS_GLOB_csvdata.py`, `CREATE_ncdf_GLOB_data_LYR.py` and `CREATE_ncdf_GLOB_data_NYA.py` serve for producing the netCDF data file "GLOB_data_5min_*.nc". You can download this file on GIT so you won't need to use these files.

1/ You need to produce the daily data from `GLOB_estimations_Faiman_LYR.py` or `GLOB_estimations_Faiman_NYA.py` with the data "GLOB_data_5min_*.nc". This can take a bit of time if one want to generate data over the whole year. It takes ca. 1h to generate data for a month.

2/ Use the newly produced data in `CALC_irradiance_multiorient.py` then `PLOT_irradiance_multiorient_polar_plot.py`.

/!\ Don't forget to change the datapath in the scripts. 


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
Arthur G - [arthurg@unis.no](mailto:arthurg@unis.no)

Project Link: [https://github.com/ArthurGarreau/GLOB](https://github.com/ArthurGarreau/GLOB)