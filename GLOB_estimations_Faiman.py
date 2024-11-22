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
os.chdir('C:/Users/arthurg/OneDrive - NTNU/Workspace/PAPER_2_Data_Analysis/GLOB_scripts')

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


# %%  Beam irradiance calculation
from itertools import combinations
import glob_functions as fct
date = '2024-04-14'

df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()

# Define the 25 variables of interest
variables = ['GHI', 'S_45', 'S_90', 'S_135', 'SW_45', 'SW_90',
             'SW_135', 'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135',
             'N_45', 'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45',
             'E_90', 'E_135', 'SE_45', 'SE_90', 'SE_135']

# Generate all possible combinations of 2 variables
combinations_2 = list(combinations(variables, 2))

# Create a dictionary to store the new columns
new_columns = {}

# Calculate new values for each combination
for var1, var2 in combinations_2:
    # Example operation: sum of the two variables
    new_columns[f"{var1}_plus_{var2}"] = np.array([df_glob_one_day[var1] , df_glob_one_day[var2]])

# Convert the dictionary to a DataFrame
new_columns_df = pd.DataFrame(new_columns, index=df_glob_one_day.index)

# Concatenate the original DataFrame with the new columns
df_combinations = pd.concat([df_glob_one_day, new_columns_df], axis=1)

# Check the resulting DataFrame
df_combinations.head()
# Drop the columns
df_combinations = df_combinations.drop(columns=df_glob_one_day.columns)

# Check the resulting DataFrame
for d in df_combinations.columns:
    print(d)


