# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 14:08:11 2025

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

kz_datafile_path = DATA_PATH.parent / "Irradiance_ncdf" / "Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"
glob_estim_path = SCRIPT_PATH / "Data" / "B_and_D_Estimations_LYR"
glob_nc_datafile = DATA_PATH / f"GLOB_data_{f}min_2023-24.nc"

output_plot_path = SCRIPT_PATH / "Figures_all" / "Beam_and_Diffuse" / "LYR"

###############################################################################

year=2024
# Define the date range
# start_date = f'{year}-04-03'
# end_date = f'{year}-04-03'
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Load the GHI data from the NetCDF file
ds = xr.open_dataset(glob_nc_datafile)
ghi_data = ds['GHI']


for date in dates:    # Generate the file paths for the current date
    file_date_str = date.strftime("%Y-%m-%d")
    glob_filename = glob_estim_path / f"estimation_beam_diffuse_{file_date_str}_{f}min.csv" 

    # Check if the GLOB file exists
    if not os.path.exists(glob_filename):
        print(f"GLOB file not found for {file_date_str}")
        continue  # Skip to the next iteration


    # Load the GLOB data
    try:
        glob_data = pd.read_csv(glob_filename, sep='\t', parse_dates=['Timestamp'], index_col='Timestamp', skiprows=5)
    except Exception as e:
        print(f"Error loading GLOB data for {file_date_str}: {e}")
        continue  # Skip to the next iteration

    # Filter GLOB data for the current date and time range
    # glob_data = glob_data[(glob_data.index.hour >= 9) & (glob_data.index.hour < 16)]

    # Extract the Beam and Diffuse components
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