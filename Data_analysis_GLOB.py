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
os.chdir('C:/Users/arthurg/OneDrive - NTNU/Workspace/PAPER_2_Data_Analysis/')

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
"""
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
"""
# %%  Beam irradiance calculation
import functions as fct
date = '2024-04-14'

df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()


# ds_kz = ds_kz.sortby('time')
# ds_kz_one_day = ds_kz.sel(time='2024-04-14')

df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')
time_series = df_glob_one_day.index

# Latitude and longitude for location (example: Longyearbyen)
latitude = ds_glob.latitude
longitude = ds_glob.longitude


solar_position = pvlib.solarposition.get_solarposition(time_series, latitude, longitude)

# Extract solar angles in radians
solar_zenith = np.radians(solar_position['zenith'].values)  # Zenith angle
solar_azimuth = np.radians(solar_position['azimuth'].values)  # Azimuth angle
hour_angle = np.radians(pvlib.solarposition.hour_angle(time_series, longitude, solar_position.equation_of_time))
sin_alpha = np.sin(np.radians(solar_position.elevation[solar_position.elevation>0]))

delta = pvlib.solarposition.declination_spencer71(time_series.day_of_year)
omega = hour_angle
phi = np.radians(latitude)
beta = np.radians(45)  # Tilt angle in radians
rho = 0.7            # Albedo (reflectivity)
gamma = {
    "East": np.radians(90),
    "South": np.radians(180),
    "West":  np.radians(270),
    "North": np.radians(0)
    }
# Directional angles for clarity
directions = ["East", "South", "West", "North"]

# Updated DataFrames with the North component
C_B_df = pd.DataFrame({
    "East": fct.C_B_lambda(phi, beta, delta, omega, gamma["East"]),
    "South": fct.C_B_lambda(phi, beta, delta, omega, gamma["South"]),
    "West": fct.C_B_lambda(phi, beta, delta, omega, gamma["West"]),
    "North": fct.C_B_lambda(phi, beta, delta, omega, gamma["North"])
}, index=time_series)

C_D_df = pd.DataFrame({
    "East": fct.C_D_lambda(phi, beta, delta, omega, gamma["East"]),
    "South": fct.C_D_lambda(phi, beta, delta, omega, gamma["South"]),
    "West": fct.C_D_lambda(phi, beta, delta, omega, gamma["West"]),
    "North": fct.C_D_lambda(phi, beta, delta, omega, gamma["North"])
}, index=time_series)

C_R_df = pd.DataFrame({
    "East": fct.C_R_lambda(beta),
    "South": fct.C_R_lambda(beta),
    "West": fct.C_R_lambda(beta),
    "North": fct.C_R_lambda(beta)
}, index=time_series)

# # %%


# # Extract coefficients for each direction
# C_B_E = C_B_df["East"]
# C_B_S = C_B_df["South"]
# C_B_W = C_B_df["West"]
# C_B_N = C_B_df["North"]

# C_D_E = C_D_df["East"]
# C_D_S = C_D_df["South"]
# C_D_W = C_D_df["West"]
# C_D_N = C_D_df["North"]

# C_R_E = C_R_df["East"]
# C_R_S = C_R_df["South"]
# C_R_W = C_R_df["West"]
# C_R_N = C_R_df["North"]
    
# denominators =pd.DataFrame( {
#     "SE" : C_D_S * C_B_E - C_B_S * C_D_E,
#     "SW" : C_D_W * C_B_S - C_B_W * C_D_S,
#     "NW" : C_D_N * C_B_W - C_B_N * C_D_W,
#     "NE" : C_D_E * C_B_N - C_B_E * C_D_N}
#     , index=time_series)

# denominators.plot()


# # %%

Irradiance =  pd.DataFrame({
    "Beam GLOB": -fct.I_B_n(df_glob_one_day, C_B_df, C_D_df, C_R_df, rho),
    "Beam NYA": df_NYA['DIR [W/m**2]'].loc[date],
    "Diffuse GLOB": fct.I_D(df_glob_one_day, phi, np.radians(90), delta, omega, gamma),
    "Diffuse NYA": df_NYA['DIF [W/m**2]'].loc[date],
    "Reflected GLOB": fct.I_R(df_glob_one_day, np.radians(90), rho),
    "GHI": df_glob_one_day['GHI']
    }).dropna()


Irradiance["Total GLOB"] =  sin_alpha * Irradiance["Beam GLOB"] + Irradiance["Diffuse GLOB"] 




ax = Irradiance.plot(figsize=(10, 6))  # Adjust the figure size as needed
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout();plt.show()
plt.ylim([0, 1200])


# plt.savefig(date + '.png', dpi=300)

# plt.plot(C_B_E, label='E')
# plt.plot(C_B_W, label='W')
# plt.plot(C_B_S, label='S')
# plt.legend()

