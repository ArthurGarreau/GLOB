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
# %% Ny Ålesund

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pvlib
import xarray as xr
from config_path import SCRIPT_PATH, DATA_PATH

# Parameters
method = 'linear'; year = 2025; f = 10 #min data frequency
pyrano_var = np.zeros(25)

############################## File Paths #####################################
B_D_estimations_datafile = SCRIPT_PATH / "Data" / "Beam_Diffuse_Estimations" / \
    f"{year}_estimation_beam_diffuse_{f}min_{method}_{len(pyrano_var)}pyrano.csv"
bsrn_datafile = DATA_PATH.parent / "Irradiance" / "NYA" / "NYA_radiation_2025-all.tab"

output_plot_path = SCRIPT_PATH / "Figures_all" / "Beam_and_Diffuse" / f'{year}'
###############################################################################
def round_up_to_hundred(x):
    return np.ceil(x / 100) * 100

# Define the date range
start_date = f'{year}-03-16'
end_date = f'{year}-03-31'
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Load the BSRN data
bsrn_data_full = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')
bsrn_data_full = bsrn_data_full.resample(f'{f}min').first()

# Load GLOB data
glob_estim_data = pd.read_csv(B_D_estimations_datafile, parse_dates=['Timestamp'], index_col='Timestamp', sep='\t', header=10)

# Load GHI data from NetCDF file
ds = xr.open_dataset(DATA_PATH / f"GLOB_data_10min_{year}.nc")
ghi_day = ds['GHI'].to_dataframe()
ds.close()

for date in dates:
    # Generate the file paths for the current date
    file_date_str = date.strftime('%Y-%m-%d')

    # Filter BSRN data for the current date and time range
    bsrn_data = bsrn_data_full.loc[file_date_str]

    # Filter GLOB data for the current date and time range
    glob_data = glob_estim_data.loc[file_date_str]

    # Extract the Beam and Diffuse components
    beam_bsrn = bsrn_data['DIR']
    diffuse_bsrn = bsrn_data['DIF']
    beam_glob = glob_data['Beam']
    diffuse_glob = glob_data['Diffuse']

    # Filter GHI data for the current date and time range
    ghi_day_current = ghi_day.loc[file_date_str]

    # Calculate beam and diffuse irradiance using the decomposition models
    time = pd.date_range(file_date_str, periods=round(60/f*24), freq=f"{f}min")
    solar_position = pvlib.solarposition.get_solarposition(time, ds.latitude.values, ds.longitude.values)
    zenith = solar_position['zenith'].values

    erbs = pvlib.irradiance.erbs(ghi_day_current['GHI'].values, zenith, time)
    beam_erbs, diffuse_erbs = erbs['dni'], erbs['dhi']

    # Calculate beam and diffuse irradiance using the DIRINT model
    perez_dni = pvlib.irradiance.dirint(ghi_day_current['GHI'].values, zenith, time)
    beam_perez, diffuse_perez = perez_dni, ghi_day_current['GHI'].values - np.cos(np.radians(zenith)) * perez_dni

    # Calculate beam and diffuse irradiance using the Orgill_and_holland model
    orgill_holland = pvlib.irradiance.orgill_hollands(ghi_day_current['GHI'].values, zenith, time)
    beam_oh, diffuse_oh = orgill_holland['dni'], orgill_holland['dhi']

    # Plot the data
    plt.rc('font', size=13)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    
    # Beam subplot (top)
    ax1.plot(beam_bsrn.index, beam_bsrn, 'k', marker='.', linewidth=2, 
             label='BSRN Ny-Ålesund (truth)')
    ax1.plot(beam_glob.index, beam_glob, 'royalblue', linewidth=2, 
             label=f'GLOB (MPA)\n{len(pyrano_var)} pyranometers {method}')
    ax1.plot(ghi_day_current.index, beam_perez,  color='lightcoral', linestyle='--', 
             label='Model Perez')
    ax1.plot(ghi_day_current.index, beam_erbs, color='firebrick', linestyle='--', 
             label='Model Erbs')
    ax1.plot(ghi_day_current.index, beam_oh,  color='maroon', linestyle='--', 
             label='Model Orgill and Hollands')
    ax1.set_title('Beam', fontsize=14)
    ax1.set_ylabel('Irradiance / $W \ m^{-2}$')
    ax1.legend()
    ax1.set_ylim((0, round_up_to_hundred(np.nanmax([beam_bsrn, beam_glob, beam_erbs, beam_oh, beam_perez]))))
    ax1.grid(True, linestyle=':')
    
    # Diffuse subplot (bottom)
    ax2.plot(diffuse_bsrn.index, diffuse_bsrn, 'k', marker='.')
    ax2.plot(diffuse_glob.index, diffuse_glob, 'royalblue', linewidth=2)
    ax2.plot(ghi_day_current.index, diffuse_perez,  color='lightcoral', linestyle='--')
    ax2.plot(ghi_day_current.index, diffuse_erbs,  color='firebrick', linestyle='--')  
    ax2.plot(ghi_day_current.index, diffuse_oh,  color='maroon', linestyle='--')
    ax2.set_title('Diffuse', fontsize=14)
    ax2.set_xlabel('Time (UTC)')
    ax2.set_ylabel('Irradiance / $W \ m^{-2}$')
    ax2.set_ylim((0, round_up_to_hundred(np.nanmax([diffuse_bsrn, diffuse_glob, diffuse_erbs, diffuse_oh, diffuse_perez]))))
    ax2.grid(True, linestyle=':')
    
    # Format x-axis
    start_time = pd.Timestamp(beam_bsrn.index[0].date())  # Today at 00:00
    end_time = start_time + pd.Timedelta(days=1)          # Tomorrow at 00:00 (24:00)
    ax1.set_xlim(start_time, end_time)
    ax2.set_xlim(start_time, end_time)
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=1))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    ax2.annotate(beam_bsrn.index[0].strftime('%Y-%m-%d'), 
             xy=(0.90, -0.1), xycoords='axes fraction')
    plt.xticks()
    plt.yticks()
    
    plt.tight_layout()
    # plt.show()
    # Uncomment the following line to save the plot
    output_plot_path.mkdir(parents=True, exist_ok=True)
    output_file = output_plot_path / f'Beam_Diffuse_{file_date_str}_{len(pyrano_var)}pyrano_{method}.png'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
# %% Longyearbyen (2023-2024)


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
import numpy as np
import pvlib
from config_path import SCRIPT_PATH, DATA_PATH

method = 'linear'; year = 2024; f = 10 #min data frequency
pyrano_var = np.zeros(5)

############################## File Paths #####################################
kz_datafile = DATA_PATH.parent / "Irradiance_ncdf" / "Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"
glob_nc_datafile = DATA_PATH / f"GLOB_data_{f}min_{year}.nc"
B_D_estimations_datafile = SCRIPT_PATH / "Data" / "Beam_Diffuse_Estimations" / \
    f"{year}_estimation_beam_diffuse_{f}min_{method}_{len(pyrano_var)}pyrano.csv"

output_plot_path = SCRIPT_PATH / "Figures_all" / "Beam_and_Diffuse" / f'{year}'
###############################################################################
def round_up_to_hundred(x):
    return np.ceil(x / 100) * 100

# Define the date range
start_date = f'{year}-05-1'
end_date = f'{year}-05-31'
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Load the GHI data from the NetCDF file
ds = xr.open_dataset(glob_nc_datafile)
ghi_data = ds['GHI'].to_pandas();
# Load GLOB data
glob_estim_data = pd.read_csv(B_D_estimations_datafile, parse_dates=['Timestamp'],
                              index_col='Timestamp', sep='\t', header=10)

for date in dates:    # Generate the file paths for the current date
    file_date_str = date.strftime("%Y-%m-%d")


    # Filter data for the current date and time range
    glob_day_data = glob_estim_data.loc[file_date_str]
    ghi_day_data = ghi_data.loc[file_date_str]
    
    # Calculate zenith
    solar_position = pvlib.solarposition.get_solarposition(glob_day_data.index,
                                                           ds.latitude.values, ds.longitude.values)
    zenith_glob = solar_position['zenith'].values
    
    # Calculate the reconstructed GHI
    glob_day_data["Constructed Global"] = np.cos(np.deg2rad(zenith_glob)) \
        * glob_day_data['Beam'] + glob_day_data['Diffuse']


    # Calculate beam and diffuse irradiance using the decomposition models
    time = pd.date_range(file_date_str, periods=round(60/f*24), freq=f"{f}min")
    solar_position = pvlib.solarposition.get_solarposition(time, ds.latitude.values, ds.longitude.values)
    zenith = solar_position['zenith'].values

    erbs = pvlib.irradiance.erbs(ghi_day_data.values, zenith, time)
    beam_erbs, diffuse_erbs = erbs['dni'], erbs['dhi']

    # Calculate beam and diffuse irradiance using the DIRINT model
    perez_dni = pvlib.irradiance.dirint(ghi_day_data.values, zenith, time)
    beam_perez, diffuse_perez = perez_dni, ghi_day_data.values - np.cos(np.radians(zenith)) * perez_dni

    # Calculate beam and diffuse irradiance using the Orgill_and_holland model
    orgill_holland = pvlib.irradiance.orgill_hollands(ghi_day_data.values, zenith, time)
    beam_oh, diffuse_oh = orgill_holland['dni'], orgill_holland['dhi']

    
    # Plot the data
    plt.rc('font', size=14)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

    # Beam subplot (top)
    ax1.plot(glob_day_data['Beam'], 'royalblue', linewidth=2, 
             label=f'GLOB (MPA)\n{len(pyrano_var)} pyranometers {method}')
    ax1.plot(beam_perez,  color='lightcoral', linestyle='--', 
             label='Model Perez')
    ax1.plot(beam_erbs, color='firebrick', linestyle='--', 
             label='Model Erbs')
    ax1.plot( beam_oh,  color='maroon', linestyle='--', 
             label='Model Orgill and Hollands')
    ax1.set_title('Beam', fontsize=14)
    ax1.set_ylabel('Irradiance / $W \ m^{-2}$')
    ax1.legend(fontsize=11, loc=8)
    ax1.set_ylim((0, round_up_to_hundred(np.nanmax([glob_day_data['Beam'], beam_erbs, beam_oh, beam_perez]))))
    ax1.grid(True, linestyle=':')
    
    # Diffuse subplot (middle)
    ax2.plot(glob_day_data['Diffuse'], 'royalblue', linewidth=2)
    ax2.plot(diffuse_perez,  color='lightcoral', linestyle='--')
    ax2.plot(diffuse_erbs,  color='firebrick', linestyle='--')  
    ax2.plot(diffuse_oh,  color='maroon', linestyle='--')
    ax2.set_title('Diffuse', fontsize=14)
    ax2.set_ylabel('Irradiance / $W \ m^{-2}$')
    ax2.set_ylim((0, round_up_to_hundred(np.nanmax([glob_day_data['Diffuse'], diffuse_erbs, diffuse_oh, diffuse_perez]))))
    ax2.grid(True, linestyle=':')
    
    # Global subplot (bottom)
    ax3.plot(glob_day_data["Constructed Global"], 'royalblue', label=f'Constructed GHI GLOB (MPA)\n{len(pyrano_var)} pyranometers {method}')
    ax3.plot(ghi_day_data, 'k', marker='.', label='Measured GHI')
    ax3.legend(fontsize=11, loc=8)
    ax3.set_title('Global', fontsize=14)
    ax3.set_xlabel('Time (UTC)')
    ax3.set_ylabel('Irradiance / $W \ m^{-2}$')
    ax3.grid(True, linestyle=':')

    # Format x-axis
    start_time = pd.Timestamp(glob_day_data.index[0].date())  # Today at 00:00
    end_time = start_time + pd.Timedelta(days=1)          # Tomorrow at 00:00 (24:00)
    ax1.set_xlim(start_time, end_time)
    ax2.set_xlim(start_time, end_time)
    ax3.set_xlim(start_time, end_time)

    ax3.xaxis.set_major_locator(mdates.HourLocator(interval=1))
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    ax3.annotate(glob_day_data.index[0].strftime('%Y-%m-%d'), 
             xy=(0.90, -0.2), xycoords='axes fraction')
    plt.xticks()
    plt.yticks()
    
    plt.tight_layout()
    plt.grid(True, linestyle=':')

    output_plot_path.mkdir(parents=True, exist_ok=True)
    output_file = output_plot_path / f'Beam_Diffuse_{file_date_str}_{len(pyrano_var)}pyrano_{method}.png'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    # plt.show()

# Close the NetCDF file
ds.close()






