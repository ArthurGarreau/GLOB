# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 17:01:25 2024

@author: arthurg

Contact: arthurg@unis.no
"""
# %% Load libraries

import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import os
os.chdir('C:/Users/arthurg/OneDrive - Universitetssenteret på Svalbard AS/Documents/UNIS_PhD/PAPER_2/PAPER_2_Data_Analysis/GLOB_scripts')

import pvlib


# %% Load data

# GLOB data NetCDF 
file_path_glob = "C:/Users/arthurg/OneDrive - NTNU/Workspace/Data/GLOB/GLOB_data.nc"
ds_glob = xr.open_dataset(file_path_glob)

# K&Z data Adventdalen
file_path_kz = "C:/Users/arthurg/OneDrive - NTNU/Workspace/Data/Irradiance_ncdf/Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"
ds_kz = xr.open_dataset(file_path_kz)

# Ny-Ålesund data 
file_path_nya = "C:/Users/arthurg/OneDrive - NTNU/Workspace/Data/Irradiance/NYA_radiation_2006-05_2024-10.csv"
df_NYA = pd.read_csv(file_path_nya,sep='\t',skiprows=23, encoding='utf-8')
df_NYA.rename(columns={'Date/Time': 'Timestamp'}, inplace=True)
df_NYA['Timestamp'] = pd.to_datetime(df_NYA['Timestamp'], format='%Y-%m-%dT%H:%M')
df_NYA = df_NYA[df_NYA['Timestamp'].dt.year >= 2023]
df_NYA.set_index('Timestamp', inplace=True)
df_NYA.index = df_NYA.index.tz_localize('UTC')


# %% Data management 

# Convert xarray datasets to DataFrames for faster operations
df_glob = ds_glob.to_dataframe().reset_index()
df_kz = ds_kz.to_dataframe()[['SWdown']].reset_index()

# Convert timestamp columns to datetime and set as index
df_glob['Timestamp'] = pd.to_datetime(df_glob['Timestamp'])
df_kz['time'] = pd.to_datetime(df_kz['time'])
df_glob.set_index('Timestamp', inplace=True)
df_kz.set_index('time', inplace=True)

# Resample GLOB data to 5-minute intervals if it’s at a higher frequency
# Only resample if required (skip if already 5-min)
if df_glob.index.freq is None or df_glob.index.freq != '5T':
    df_glob = df_glob.resample('5min').mean()

# Merge DataFrames on timestamps (inner join for synchronization)
merged_df = pd.merge(df_glob, df_kz, left_index=True, right_index=True, how='inner').dropna()
merged_df = merged_df.resample('1D').mean().dropna()



# %% Plot regression - Figure 1


# Perform linear regression between 'GHI' and 'SWdown'
slope, intercept, r_value, p_value, std_err = linregress(merged_df['GHI'], merged_df['SWdown'])

# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(merged_df['GHI'], merged_df['SWdown'], label='Data points', color='k', alpha=0.5)
plt.plot(merged_df['GHI'], intercept + slope * merged_df['GHI'], 'r', label=f'Linear fit: y={slope:.2f}x + {intercept:.2f} | R² = {r_value**2:.3f}')
# Add the f(x) = x line
min_val = min(merged_df['GHI'].min(), merged_df['SWdown'].min())
max_val = max(merged_df['GHI'].max(), merged_df['SWdown'].max())
plt.plot([min_val, max_val], [min_val, max_val], 'k--', label='y = x')

plt.xlabel('GHI GLOB (W/m²)')
plt.ylabel('GHI Kipp & Zonen CNR1 (W/m²)')
plt.legend()
plt.grid(True)
plt.savefig('Fig 1 Linear regression daily GLOB K&Z.png', dpi=300)
plt.close()


