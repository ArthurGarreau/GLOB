# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 12:20:37 2025

@author: arthurg
"""

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from config_path import ALL_PLOT_PATH, GLOB_DATA_PATH

# Define the file path
file_path = GLOB_DATA_PATH / "GLOB_data_10min_2023-24.nc"


# Define months to analyze
months = [4, 5, 6]  # April, May, June
month_names = ['April', 'May', 'June']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Colors for each month
e_markers = ['o', 's', '^']  # Markers for E_45: circle, square, triangle
w_markers = ['o', 's', '^']  # Markers for W_45: diamond, pentagon, star

# Load the data
print("Loading data...")
ds = xr.open_dataset(file_path)
df = ds.to_dataframe()

# Convert index to timezone-aware datetime if needed
if df.index.tz is None:
    df.index = df.index.tz_localize('UTC')
df.index = df.index.tz_convert('Etc/GMT-1')  # CET without DST


# Filter for the months of interest
print("Filtering data for April, May, June...")
df_filtered = df[df.index.month.isin(months)]

# Create a figure with larger size to accommodate bigger legend
plt.figure(figsize=(14, 7))
plt.rcParams['font.size'] = 12

# Process and plot each month
for i, month in enumerate(months):
    # Filter data for the current month
    month_data = df_filtered[df_filtered.index.month == month]

    # Group by hour and calculate mean and standard deviation
    e_45_hourly = month_data.groupby(pd.cut(month_data.index.hour + month_data.index.minute/60,
                                           bins=24, labels=False))['E_45'].agg(['mean', 'std'])
    w_45_hourly = month_data.groupby(pd.cut(month_data.index.hour + month_data.index.minute/60,
                                            bins=24, labels=False))['W_45'].agg(['mean', 'std'])
        
    # Create time axis for plotting (centered on each hour)
    hours_centered = np.arange(0.5, 24, 1)

    # Plot E_45 with solid line and specific marker
    plt.plot(hours_centered, e_45_hourly['mean'], color=colors[i], linestyle='-',
             label=f'E_45 ({month_names[i]})', linewidth=2,
             marker=e_markers[i], markersize=8)

    # Add shaded area for standard deviation
    plt.fill_between(hours_centered,
                     e_45_hourly['mean'] - e_45_hourly['std'],
                     e_45_hourly['mean'] + e_45_hourly['std'],
                     color=colors[i], alpha=0.15)

    # Plot W_45 with dashed line and specific marker
    plt.plot(hours_centered, w_45_hourly['mean'], color=colors[i], linestyle='--',
             label=f'W_45 ({month_names[i]})', linewidth=2,
             marker=w_markers[i], markersize=9)

    # Add shaded area for standard deviation
    plt.fill_between(hours_centered,
                     w_45_hourly['mean'] - w_45_hourly['std'],
                     w_45_hourly['mean'] + w_45_hourly['std'],
                     color=colors[i], alpha=0.1)



# Add plot elements
# Add vertical line at 11:00
plt.axvline(x=12, color='k', linestyle=':', linewidth=2)
plt.title('Average Daily Evolution of E_45 and W_45 Irradiance (April-June 2023-2024)', fontsize=16, pad=20)
plt.xlabel('Hour of Day (CET)', fontsize=14)
plt.ylabel('Irradiance [$W \ m^{-2}$]', fontsize=14)
plt.ylim(0, None)  # Let matplotlib determine the upper y-limit

# Set x-ticks to show every 2 hours
plt.xticks(np.arange(0, 24, 1))

# Add grid
plt.grid(True, linestyle=':', alpha=0.7)

# Create legend with larger font size and better positioning
plt.legend(fontsize=12, ncol=3, loc="upper left", framealpha=1, edgecolor='black')

# Adjust layout to prevent label cutoff and make room for legend
plt.tight_layout(rect=[0, 0.08, 1, 0.95])

# Save the plot
plt.savefig(ALL_PLOT_PATH / "E_W_45_daily_evolution.png", dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

# Close the dataset
ds.close()
print("Done!")
