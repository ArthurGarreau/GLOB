# -*- coding: utf-8 -*-
"""
Beam and Diffuse Components Plotting Script
============================================

This script plots the Beam and Diffuse components from Ny-Ålesund BSRN measurements
and Ny-Ålesund GLOB instrument data on the same plot for each available day between 9:00 and 16:00 UTC.
It also includes the GHI from the NetCDF file and the reconstructed GHI.

Key Features:
-------------
- Loads BSRN and GLOB data from specified file paths.
- Parses datetime formats correctly.
- Filters data by date and time range (9:00 to 16:00 UTC).
- Plots Beam, Diffuse, GHI, and reconstructed GHI components on the same graph.
- Saves the plots as PNG files with a detailed header including metadata.

Dependencies:
-------------
- numpy
- pandas
- matplotlib
- pathlib
- datetime
- os
- xarray

Created on Mon Apr 28 16:23:23 2025

@author: arthurg
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime, timedelta
import os
import xarray as xr
import numpy as np

############################## File Paths #####################################
bsrn_datafile = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\Irradiance\NYA") / "NYA_radiation_2025-03.tab"
glob_estim_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB\B_and_D_Estimations_NYA")
glob_nc_datafile = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB") / "GLOB_data_5min_2025.nc"

output_plot_path = Path(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\Fig\Estim_Faiman\nya")

###############################################################################

# Define the date range
start_date = datetime(2025, 3, 16)
end_date = datetime(2025, 4, 26)

# Load the GHI data from the NetCDF file
ds = xr.open_dataset(glob_nc_datafile)
ghi_data = ds['GHI']

# Loop through each day
current_date = start_date
while current_date <= end_date:
    # Generate the file paths for the current date
    file_date_str = current_date.strftime("%Y-%m-%d")
    glob_filename = glob_estim_path / f"best_estimations_{file_date_str}_ERBS.csv"

    # Check if the GLOB file exists
    if not os.path.exists(glob_filename):
        print(f"GLOB file not found for {file_date_str}")
        current_date += timedelta(days=1)
        continue

    # Load the BSRN data
    try:
        bsrn_data = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')
    except Exception as e:
        print(f"Error loading BSRN data: {e}")
        current_date += timedelta(days=1)
        continue

    # Filter BSRN data for the current date and time range
    bsrn_data = bsrn_data[(bsrn_data.index.normalize() == pd.Timestamp(file_date_str))]
    # bsrn_data = bsrn_data[(bsrn_data.index.hour >= 9) & (bsrn_data.index.hour < 16)]

    # Load the GLOB data
    try:
        glob_data = pd.read_csv(glob_filename, sep='\t', parse_dates=['Timestamp'], index_col='Timestamp', skiprows=5)
    except Exception as e:
        print(f"Error loading GLOB data for {file_date_str}: {e}")
        current_date += timedelta(days=1)
        continue

    # Filter GLOB data for the current date and time range
    # glob_data = glob_data[(glob_data.index.hour >= 9) & (glob_data.index.hour < 16)]

    # Extract the Beam and Diffuse components
    beam_bsrn = bsrn_data['DIR [W/m**2]']
    diffuse_bsrn = bsrn_data['DIF [W/m**2]']
    beam_glob = glob_data['Beam']
    diffuse_glob = glob_data['Diffuse']

    # Calculate the reconstructed GHI
    glob_data["Constructed Global"] = np.cos(np.deg2rad(glob_data['Solar zenith'])) * glob_data['Beam'] + glob_data['Diffuse']

    # Filter GHI data for the current date and time range
    ghi_day = ghi_data.sel(Timestamp=f"{file_date_str}")

    # Plot the data
    plt.figure(figsize=(12, 8))
    plt.plot(beam_bsrn.index, beam_bsrn, 'k', label='Beam BSRN')
    plt.plot(diffuse_bsrn.index, diffuse_bsrn, 'k:', label='Diffuse BSRN')
    plt.plot(beam_glob.index, beam_glob, 'r', alpha=0.7, label='Beam GLOB')
    plt.plot(diffuse_glob.index, diffuse_glob, 'r:', alpha=0.7, label='Diffuse GLOB')
    plt.plot(ghi_day.Timestamp, ghi_day, 'g', label='GHI pyr GLOB')
    plt.plot(glob_data.index, glob_data["Constructed Global"], 'g--', label='Constructed GHI')

    plt.xlabel('Time (UTC)')
    plt.ylabel('Irradiance (W m-2)')
    plt.title(f'Beam and Diffuse Components - Ny-Ålesund - {file_date_str}')
    plt.legend()
    plt.grid(True, linestyle=':')

    # Save the plot
    output_plot_path.mkdir(parents=True, exist_ok=True)
    output_file = output_plot_path / f"beam_diffuse_plot_{file_date_str}.png"
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

    print(f"Plot for {file_date_str} saved to {output_file}")

    # Move to the next day
    current_date += timedelta(days=1)

# Close the NetCDF file
ds.close()
