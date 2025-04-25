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

### `CALC_irradiance_multiorient.py`
This script calculates the irradiance for multiple orientations based on the data collected by the GLOB instrument.

### `PLOT_irradiance_multiorient_polar_plot.py`
This script generates polar plots to visualize the irradiance data from different orientations.

### `Data_analysis_GLOB.py`
This script performs various data analysis tasks on the GLOB dataset, including cleaning, processing, and summarizing the data.

### `GLOB_Polar_Heatmap.py`
This script creates polar heatmaps to visualize the distribution of solar irradiance across different orientations and times.

### `GLOB_estimations_Faiman_LYR.py`
This script calculates beam and diffuse irradiance using the Faiman et al. (1992) method. It processes GLOB and K&Z data, estimates irradiance components, and generates plots.

### `GLOB_estimations_Faiman_NYA.py`
This script calculates beam and diffuse irradiance using the Faiman et al. (1992) method. It processes GLOB and K&Z data, estimates irradiance components, and generates plots.

### `glob_functions_Faiman.py`
This script contains functions for calculating solar angles, irradiance components, and geometry coefficients. It uses the Faiman et al. (1992) method for beam and diffuse irradiance estimation.

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

To run the scripts, use the following commands:

```bash
python CALC_irradiance_multiorient.py
python PLOT_irradiance_multiorient_polar_plot.py
python Data_analysis_GLOB.py
python GLOB_Polar_Heatmap.py
python GLOB_estimations_Faiman_LYR.py
python GLOB_estimations_Faiman_NYA.py
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
Arthur G - [arthurg@unis.no](mailto:arthurg@unis.no)

Project Link: [https://github.com/ArthurGarreau/GLOB](https://github.com/ArthurGarreau/GLOB)