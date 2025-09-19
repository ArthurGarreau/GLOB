# GLOB: Multiface Solar Irradiance Analysis Toolkit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**GLOB** is a Python-based toolkit for analyzing solar irradiance data from a 26-pyranometer array (Apogee SP-110), designed for Arctic environments. This project is part of my PhD research at The University Centre in Svalbard, funded by The Arctic Field Grant (Grant nr. xxx).

---

## Table of Contents
- [Features](#features)
- [Pyranometer Orientations](#pyranometer-orientations)
- [Scripts](#scripts)
- [Data Workflow](#data-workflow)
- [Installation](#installation)
- [Usage](#usage)
- [Outputs](#outputs)
- [License](#license)
- [Contact](#contact)

---

## Features
- **Comprehensive Irradiance Analysis**: Measures global tilted irradiance (GTI), beam, diffuse, and reflected components from 26 orientations.
- **Advanced Estimation Methods**: Implements the Faiman et al. (1987) method for beam/diffuse separation, with support for linear/nonlinear least squares and Monte Carlo simulations.
- **Data Processing Pipeline**: Converts raw CSV data to NetCDF, calculates solar angles, and estimates irradiance components.
- **Visualization Tools**: Generates polar plots, time series, and error propagation visualizations.
- **Parallel Processing**: Optimized for efficient computation using `joblib`.
- **Error Propagation**: Includes scripts for uncertainty quantification and best-combination selection.

---

## Pyranometer Orientations

The following table describes the **variable names in the NetCDF data files** for GLOB measurements, along with their azimuth and tilt angles.

| Original Name | NetCDF Variable | Azimuth (°) | Tilt (°) | Comments                |
|---------------|-----------------|-------------|----------|-------------------------|
| SIR(1)        | GHI             | 0           | 0        | Horizontal plane (GHI)  |
| SIR(2)        | S_45            | 180         | 45       | South-facing            |
| SIR(3)        | S_90            | 180         | 90       | South-facing            |
| SIR(4)        | S_135           | 180         | 135      | South-facing            |
| SIR(5)        | SW_45           | 225         | 45       | Southwest-facing        |
| SIR(6)        | SW_90           | 225         | 90       | Southwest-facing        |
| SIR(7)        | SW_135          | 225         | 135      | Southwest-facing        |
| SIR(8)        | W_45            | 270         | 45       | West-facing             |
| SIR(9)        | W_90            | 270         | 90       | West-facing             |
| SIR(10)       | W_135           | 270         | 135      | West-facing             |
| SIR(11)       | NW_45           | 315         | 45       | Northwest-facing        |
| SIR(12)       | NW_90           | 315         | 90       | Northwest-facing        |
| SIR(13)       | NW_135          | 315         | 135      | Northwest-facing        |
| SIR(14)       | N_45            | 0           | 45       | North-facing            |
| SIR(15)       | N_90            | 0           | 90       | North-facing            |
| SIR(16)       | N_135           | 0           | 135      | North-facing            |
| SIR(17)       | NE_45           | 45          | 45       | Northeast-facing        |
| SIR(18)       | NE_90           | 45          | 90       | Northeast-facing        |
| SIR(19)       | NE_135          | 45          | 135      | Northeast-facing        |
| SIR(20)       | E_45            | 90          | 45       | East-facing             |
| SIR(21)       | E_90            | 90          | 90       | East-facing             |
| SIR(22)       | E_135           | 90          | 135      | East-facing             |
| SIR(23)       | SE_45           | 135         | 45       | Southeast-facing        |
| SIR(24)       | SE_90           | 135         | 90       | Southeast-facing        |
| SIR25         | SE_135          | 135         | 135      | Southeast-facing        |
| Kglob         | GHI_ground      | 0           | 180      | Reflected GHI           |

> **Note:** The NetCDF data files use the naming convention shown in the **NetCDF Variable** column. The original column names (e.g., `SIR(1)`, `SIR(2)`, etc.) are renamed during preprocessing for clarity and consistency.

---

## Scripts

### **0. Path Scripts
- **`config_path.py`**: Central configuration file containing all file paths. Must be modified to match your system. Running this script creates all necessary project folders:
  ```bash
  python config_path.py

### **1. Data Preprocessing**
- **`DATA_GLOB_csvdata.py`**
  Processes raw CSV data, standardizes naming, converts timestamps to UTC, and outputs cleaned data (`GLOB_data_30sec_*.dat`).

- **`DATA_create_ncdf_GLOB_LYR.py` / `DATA_create_ncdf_GLOB_NYA.py`**
  Converts CSV data to NetCDF (`GLOB_data_5min_*.nc`), calculates solar angles, and computes surface albedo.

### **2. Irradiance Estimation**
- **`CALC_GTI.py`**
  Estimates **Global Tilted Irradiance (GTI)** for all pyranometer orientations using the Faiman et al. (1987) method.

- **`CALC_beam_diffuse.py`**
  Estimates beam and diffuse irradiance using the Faiman et al. (1987) method (linear/nonlinear least squares).

- **`CALC_beam_diffuse_3pyrano.py`**
  Estimates irradiance using only 3 pyranometers (for minimal configurations).

- **`CALC_beam_diffuse_bestcombination.py`**
  Finds the optimal pyranometer combination for least-squares estimation via parallel processing.

- **`CALC_beam_diffuse_error_propagation.py`**
  Quantifies uncertainty in beam/diffuse estimates using Monte Carlo simulations.

### **3. Metrics and Validation**
- **`CALC_metrics_beam_diffuse.py`**
  Computes performance metrics (RMSE, MBE, R²) for beam/diffuse estimates against reference data.

- **`CALC_metrics_GTI.py`**
  Validates GTI estimates using statistical metrics.

### **4. Visualization**
- **`PLOT_PAPER2.py`**
  Generates publication-ready plots (polar heatmaps, time series, error distributions).

### **5. Core Functions**
- **`glob_functions_calculation.py`**
  Provides reusable functions for solar angle calculations, irradiance transposition coefficients, and least-squares optimization.

---

## Data Workflow

### 1.1 Use Preprocessed Data (Recommended)
  Preprocessed NetCDF data files are available for download from the **Arctic Data Centre** (link will be added soon). If you use these files, you can **skip the data preprocessing step** and directly proceed to irradiance estimation.

### 1.2 Process Raw Data
   Run `DATA_GLOB_csvdata.py` → `DATA_create_ncdf_*.py` to generate NetCDF files.
   *(Note: This step is only necessary if you have raw CSV data.)*

### 2. Estimation
   Use `CALC_GTI.py`, `CALC_beam_diffuse*.py` to estimate GTI, beam, and diffuse irradiance from NetCDF data.
   *(~1 hour/month for full-year data.)*

### 3. Error Analysis
   Run `CALC_beam_diffuse_error_propagation.py` for uncertainty quantification.

### 4. Visualization
   Use `PLOT_PAPER2.py` to visualize results.

---

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/ArthurGarreau/GLOB.git
   cd GLOB

2. **Adapt the paths**
  Modify the path in the script "config_path.py" and run it to create the necessary folders.