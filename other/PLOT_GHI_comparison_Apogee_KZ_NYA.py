# -*- coding: utf-8 -*-
"""
Comparison of GLOB and BSRN Irradiance Data with Statistical Analysis
=====================================================================
This script compares GLOB pyranometer data with BSRN reference data for Ny-Ålesund, Svalbard.
It performs statistical analysis and visualizes the relationship between GHI, calculated GHI, and albedo from both datasets.
Key Features:
-------------
- Loads GLOB and BSRN irradiance data.
- Calculates solar position and derived quantities (e.g., GHI_calc, albedo).
- Plots KDE-based density scatter plots for GHI, calculated GHI, and albedo comparisons.
- Computes linear regression and standard deviation for each comparison.

Dependencies:
-------------
- xarray
- pandas
- numpy
- matplotlib
- scipy.stats
- Custom module: glob_functions_calculation

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: August 22, 2024
"""

# %% Load Libraries
import xarray as xr
import pandas as pd
from config_path import PLOT_PATH, GLOB_DATA_PATH, NYA_DATA_PATH
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress, gaussian_kde
import matplotlib.colors as colors

############################## File Paths #####################################

glob_datafile = GLOB_DATA_PATH / "GLOB_data_10min_2025.nc"
bsrn_datafile = NYA_DATA_PATH / "NYA_radiation_2025-all.tab"

output_file = PLOT_PATH / "GHI_Albedo_model_regressions.png"

###############################################################################
# Load GLOB data as a pandas DataFrame and localize timezone to UTC
df_glob = xr.open_dataset(glob_datafile).to_dataframe()
df_glob.index = df_glob.index.tz_localize('UTC')

# Load the BSRN data, skip metadata rows, parse dates, and localize timezone to UTC
df_bsrn = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')
df_bsrn = df_bsrn.resample('10min').first()
df_bsrn.index = df_bsrn.index.tz_localize('UTC')

# Align GLOB GHI and BSRN SWD data for comparison
df_merged = pd.merge(
    df_glob[['GHI']],
    df_bsrn[['SWD']],
    left_index=True,
    right_index=True,
    how='inner'
).dropna()

# Calculate albedo from BSRN data (SWU/SWD) and filter for midday values (10 AM to 12 PM)
albedo_bsrn = df_bsrn['SWU'] / df_bsrn['SWD']
filtered_albedo = albedo_bsrn.where((albedo_bsrn.index.hour >= 10) & (albedo_bsrn.index.hour < 12))

# Compute daily mean albedo and apply to all timestamps of the day
daily_mean_albedo = filtered_albedo.resample('1D').mean()
daily_mean_albedo[daily_mean_albedo < 0] = np.nan
new_albedo_bsrn = filtered_albedo * np.nan
for day in daily_mean_albedo.index:
    mask = new_albedo_bsrn.index.date == day.date()
    new_albedo_bsrn.loc[mask] = daily_mean_albedo.loc[day]

# Align GLOB and BSRN albedo data for comparison
df_merged_albedo = pd.merge(
    df_glob['albedo'].to_frame(name='albedo_glob'),
    filtered_albedo.to_frame(name='albedo_bsrn'),
    left_index=True,
    right_index=True,
    how='inner'
).dropna()
df_merged_albedo= df_merged_albedo.loc['2025-03-20':'2025-04-20']

# Function to plot KDE-based density scatter plots
def plot_kde_scatter(x, y, ax, **kwargs):
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    # Rescale z to the range [0, 1]
    z_normalized = (z - z.min()) / (z.max() - z.min())
    sc = ax.scatter(x, y, c=z_normalized, cmap='viridis', marker='.', **kwargs)
    return sc

# Create a figure with 2 subplots for GHI and albedo comparisons
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Subplot 1: SWD vs GHI (KDE)
sc = plot_kde_scatter(df_merged['GHI'], df_merged['SWD'], axes[0])
slope, intercept, r_value, p_value, std_err = linregress(df_merged['GHI'], df_merged['SWD'])
line = slope * df_merged['GHI'] + intercept
r_squared = r_value**2
std_diff = np.nanstd(df_merged['SWD'] - df_merged['GHI']) / np.nanmean(df_merged['SWD']) * 100
axes[0].plot(df_merged['GHI'], line, color='red', label=f'Regression: y = {slope:.2f}x + {intercept:.2f} | $R^2$ = {r_squared:.2f}')
axes[0].plot(df_merged['GHI'], df_merged['GHI'], color='k', linestyle='--', label='1:1 Line')
axes[0].set_xlabel('GLOB GHI [$W \ m^{-2}$]')
axes[0].set_ylabel('BSRN GHI [$W \ m^{-2}$]')
axes[0].set_title('GHI BSRN Vs GLOB')
axes[0].text(0.05, 0.75, f'Std Dev = {std_diff:.2f} %', transform=axes[0].transAxes, verticalalignment='top')
axes[0].legend()
axes[0].grid(linestyle=':')
fig.colorbar(sc, ax=axes[0], label='Density')

# Subplot 2: Albedo (KDE)
sc = plot_kde_scatter(df_merged_albedo['albedo_glob'], df_merged_albedo['albedo_bsrn'], axes[1])
slope, intercept, r_value, p_value, std_err = linregress(df_merged_albedo['albedo_glob'], df_merged_albedo['albedo_bsrn'])
line = slope * df_merged_albedo['albedo_glob'] + intercept
r_squared = r_value**2
std_diff = np.std(df_merged_albedo['albedo_bsrn'] - df_merged_albedo['albedo_glob'])
axes[1].plot(df_merged_albedo['albedo_glob'], line, color='red', label=f'Regression: y = {slope:.2f}x + {intercept:.2f} | $R^2$ = {r_squared:.2f}')
axes[1].plot([0,1], [0,1], color='k', linestyle='--', label='1:1 Line')
axes[1].set_xlim([0,1]); axes[1].set_ylim([0,1])
axes[1].set_xlabel('GLOB Albedo')
axes[1].set_ylabel('BSRN Albedo')
axes[1].set_title('Albedo')
axes[1].text(0.25, 0.75, f'Std Dev = {round(std_diff*100)}%', transform=axes[1].transAxes, verticalalignment='top')
axes[1].legend(loc="upper left")
axes[1].grid(linestyle=':')

fig.colorbar(sc, ax=axes[1], label='Density')
plt.tight_layout()
# plt.savefig(output_file, bbox_inches='tight')
# plt.close()


