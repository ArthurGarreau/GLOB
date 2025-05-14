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
import matplotlib.dates as mdates
import os
import xarray as xr
import numpy as np
import pvlib
from config import SCRIPT_PATH, DATA_PATH

f = 10 #min data frequency

############################## File Paths #####################################
bsrn_datafile = DATA_PATH.parent / "Irradiance" / "NYA" / "NYA_radiation_2025-03.tab"
glob_estim_path = SCRIPT_PATH / 'Data' / "B_and_D_Estimations_NYA"
# glob_estim_path = DATA_PATH / "B_and_D_Estimations_NYA"

glob_nc_datafile = DATA_PATH / f"GLOB_data_{f}min_2025.nc"

output_plot_path = SCRIPT_PATH / "Figures_all" / "Beam_and_Diffuse" / "NYA"

###############################################################################
year=2025
# Define the date range
# start_date = f'{year}-04-03'
# end_date = f'{year}-04-03'
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Load the GHI data from the NetCDF file
ds = xr.open_dataset(glob_nc_datafile)
ghi_data = ds['GHI']

# Load the BSRN data
bsrn_data_full = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')
bsrn_data_full = bsrn_data_full.resample(f'{f}min').first()


for date in dates:    # Generate the file paths for the current date
    file_date_str = date.strftime("%Y-%m-%d")
    glob_filename = glob_estim_path / f"estimation_beam_diffuse_{file_date_str}_{f}min.csv" 

    # Check if the GLOB file exists
    if not os.path.exists(glob_filename):
        print(f"GLOB file not found for {file_date_str}")
        continue  # Skip to the next iteration


    # Filter BSRN data for the current date and time range
    bsrn_data = bsrn_data_full[(bsrn_data_full.index.normalize() == pd.Timestamp(file_date_str))]
    # bsrn_data = bsrn_data[(bsrn_data.index.hour >= 9) & (bsrn_data.index.hour < 16)]

    # Load the GLOB data
    try:
        glob_data = pd.read_csv(glob_filename, sep='\t', parse_dates=['Timestamp'], index_col='Timestamp', skiprows=5)
    except Exception as e:
        print(f"Error loading GLOB data for {file_date_str}: {e}")
        continue  # Skip to the next iteration

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
    
    # Calculate beam and diffuse irradiance using the Erbs model
    time = pd.date_range(date, periods=round(60/f*24), freq=f"{f}min")

    solar_position = pvlib.solarposition.get_solarposition(time, ds.latitude.values, ds.longitude.values)
    zenith = solar_position['zenith'].values
    ERBS  = pvlib.irradiance.erbs(ghi_day.values, zenith, time)
    beam_erbs, diffuse_erbs = ERBS['dni'], ERBS['dhi'] 
    
    # Plot the data
    plt.figure(figsize=(12, 8))
    plt.plot(beam_bsrn.index, beam_bsrn, 'k', label='Beam BSRN')
    plt.plot(diffuse_bsrn.index, diffuse_bsrn, 'k:', label='Diffuse BSRN')
    plt.plot(beam_glob.index, beam_glob, 'r', alpha=0.7, label='Beam GLOB')
    plt.plot(diffuse_glob.index, diffuse_glob, 'r:', alpha=0.7, label='Diffuse GLOB')
    plt.plot(ghi_day.Timestamp, beam_erbs, 'b', label='Beam ERBS')
    plt.plot(ghi_day.Timestamp, diffuse_erbs, 'b:', label='Diffuse ERBS')
    # plt.plot(ghi_day.Timestamp, ghi_day, 'g', label='GHI pyr GLOB')
    # plt.plot(glob_data.index, glob_data["Constructed Global"], 'g--', label='Constructed GHI')

    plt.xlabel('Time (UTC)')
    plt.ylabel('Irradiance (W m-2)')
    plt.title(f'Beam and Diffuse Components - Ny-Ålesund - {file_date_str}')
    plt.legend()
    # Set hourly ticks
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=1))  # Tick every hour
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))  # Format as HH:MM
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim((0,800))
    plt.grid(True, linestyle=':')

    output_plot_path.mkdir(parents=True, exist_ok=True)
    output_file = output_plot_path / f"beam_diffuse_plot_{file_date_str}.png"
    # plt.savefig(output_file, bbox_inches='tight')
    plt.show()

# Close the NetCDF file
ds.close()
