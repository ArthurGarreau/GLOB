# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:35:35 2024

@author: arthurg

Contact: arthurg@unis.no
"""

# Method inspired by Faiman et al. (1992)
# Faiman, D., D. Feuermann, P. Ibbetson, and A. Zemel, 1992: A multipyranometer instrument for obtaining the solar beam and diffuse components, and the irradiance on inclined planes. Solar Energy, 48, 253–259, https://doi.org/10.1016/0038-092X(92)90099-V.


# %% Load libraries

import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import os
from pathlib import Path
import pvlib
import sys
sys.path.append(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\GLOB_scripts")


# %% Load data
data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data")

# GLOB data NetCDF 
ds_glob = xr.open_dataset(data_path / r"GLOB\GLOB_data_angles_5min.nc")
# K&Z data Adventdalen
ds_kz = xr.open_dataset(data_path / r"Irradiance_ncdf\Adventdalen_global_horizontal_irradiances_LW_SW_all.nc")

# Ny-Ålesund data 
file_path_nya = data_path / r"Irradiance\NYA_radiation_2006-05_2024-10.csv"
df_NYA = pd.read_csv(file_path_nya,sep='\t',skiprows=23, encoding='utf-8')
df_NYA.rename(columns={'Date/Time': 'Timestamp'}, inplace=True)
df_NYA['Timestamp'] = pd.to_datetime(df_NYA['Timestamp'], format='%Y-%m-%dT%H:%M')
df_NYA = df_NYA[df_NYA['Timestamp'].dt.year >= 2023]
df_NYA.set_index('Timestamp', inplace=True)
df_NYA.index = df_NYA.index.tz_localize('UTC')


# %%  Beam irradiance calculation
from itertools import combinations
import glob_functions_Faiman as fct

dates = pd.date_range(start="2024-05-07", end="2024-05-15", freq="D", tz="UTC")


for date in dates: 
    date=date.strftime('%Y-%m-%d')
    
    df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()
    df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')
    
    
    # Define the 25 variables of interest
    variables = ['GHI', 'S_45', 'S_90', 'S_135', 'SW_45', 'SW_90',
                 'SW_135', 'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135',
                 'N_45', 'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45',
                 'E_90', 'E_135', 'SE_45', 'SE_90', 'SE_135']
    
    # Generate all possible combinations of 2 variables
    combinations_2 = list(combinations(variables, 2))
    
    orientations_dict = fct.create_variable_dict(variables)
    
    # Prepare to store the results
    results = []
    
    # Loop over each minute in the day
    timestamps = pd.date_range(start=f"{date} 07:00:00", end=f"{date} 17:00:00", freq='10min', tz='UTC')
    for timestamp in timestamps:
        timestamp = pd.DatetimeIndex([timestamp])  # Ensure timestamp is in the desired format
        glob_value = df_glob_one_day.loc[timestamp]
        NYA_value = df_NYA.loc[timestamp]
    
        # Find the best combination and corresponding estimates
        best_combination, best_X_estim = fct.find_best_estimation(
            combinations_2, glob_value, NYA_value, orientations_dict, ds_glob
        )
    
        if best_X_estim is not None:
            # Append the results for this timestamp
            results.append({
                'Timestamp': timestamp[0],
                'Beam': best_X_estim[0][0],  # Beam value
                'Diffuse': best_X_estim[1][0],  # Diffuse value
                'Best_Combination': best_combination,
                'NYA_Beam': NYA_value['DIR [W/m**2]'].iloc[0],
                'NYA_Diffuse': NYA_value['DIF [W/m**2]'].iloc[0]
            })
        else:
            # Handle cases where no valid estimation is found
            results.append({
                'Timestamp': timestamp[0],
                'Beam': np.nan,
                'Diffuse': np.nan,
                'Best_Combination': np.nan,
                'NYA_Beam': NYA_value['DIR [W/m**2]'].iloc[0] if not NYA_value['DIR [W/m**2]'].isna().iloc[0] else np.nan,
                'NYA_Diffuse': NYA_value['DIF [W/m**2]'].iloc[0] if not NYA_value['DIF [W/m**2]'].isna().iloc[0] else np.nan
            })
        
        
    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    
    results_df.to_csv(data_path / "GLOB" / f"best_estimations_{date}.csv", index=False)

# # Accessing values
# print(orientations_dict['S_45']['Azimuth'])  # Output: 180
# print(orientations_dict['NE_90']['Beta'])  


# %% Plot Beam, Diffuse, NYA_Beam, and NYA_Diffuse over time

import matplotlib.dates as mdates

dates = pd.date_range(start="2024-05-07", end="2024-05-15", freq="D", tz="UTC")

for date in dates:
        
    date=date.strftime('%Y-%m-%d')
    results_df = pd.read_csv(data_path / "GLOB" / f'best_estimations_{date}.csv')  # Adjust path as needed
    
    results_df['Timestamp'] = pd.to_datetime(results_df['Timestamp'])
    
    plt.figure(figsize=(9, 6))
    
    # Plot Beam
    plt.plot(results_df['Timestamp'], results_df['Beam'], label='Longyearbyen Beam', color='orange')
    # Plot Diffuse
    plt.plot(results_df['Timestamp'], results_df['Diffuse'], label='Longyearbyen Diffuse', color='blue')
    # Plot NYA_Beam
    plt.plot(results_df['Timestamp'], results_df['NYA_Beam'], label='Ny-Ålesund Beam', color='orange', linestyle='--')
    # Plot NYA_Diffuse
    plt.plot(results_df['Timestamp'], results_df['NYA_Diffuse'], label='Ny-Ålesund Diffuse', color='blue', linestyle='--')
    
    # Customize the plot
    plt.title(f'Beam and Diffuse Components on {date}', fontsize=16)
    plt.xlabel('Time (UTC)', fontsize=14)
    plt.ylabel('Irradiance (W/m²)', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True)
    
    # Set hourly ticks
    plt.gca().xaxis.set_major_locator(mdates.HourLocator())  # Tick every hour
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))  # Format as HH:MM
    plt.xticks(fontsize=12, rotation=45)
    plt.yticks(fontsize=12)
    
    # Show the plot
    plt.tight_layout()
    
    save_path = "C:/Users/arthurg/OneDrive - Universitetssenteret på Svalbard AS/Documents/UNIS_PhD/PAPER_2/PAPER_2_Data_Analysis/Fig/"
    plt.savefig(save_path + f"Estim_Faiman_{date}")
    plt.close()
    
    
