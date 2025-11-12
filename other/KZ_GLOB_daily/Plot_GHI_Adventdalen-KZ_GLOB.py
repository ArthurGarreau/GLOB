# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 15:13:17 2025

@author: arthurg
"""

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from calendar import month_abbr
import numpy as np
import pandas as pd

# Paths to your datasets
file_kz = r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\Irradiance_ncdf\Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"
file_glob = r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB\GLOB_data_10min_2023-24.nc"

# Open the datasets
ds_kz = xr.open_dataset(file_kz)
ds_glob = xr.open_dataset(file_glob)

# Extract GHI from each dataset and convert to DataFrame
df_kz = ds_kz['SWdown'].to_dataframe().reset_index()
df_glob = ds_glob['GHI'].to_dataframe().reset_index()

# Rename columns for clarity
df_kz = df_kz.rename(columns={'time': 'timestamp', 'SWdown': 'GHI'})
df_glob = df_glob.rename(columns={'Timestamp': 'timestamp'})

# Set timestamp as index for resampling
df_kz = df_kz.set_index('timestamp')
df_glob = df_glob.set_index('timestamp')

# Resample (if needed)
df_glob.index = df_glob.index - pd.Timedelta(minutes=30)
df_glob = df_glob.resample('1h', label='right', closed='right').mean()


df_kz.index = df_kz.index - pd.Timedelta(minutes=30)
df_kz = df_kz.resample('1h', label='right', closed='right').mean()

df_kz = df_kz.loc['2023-01-01':'2024-12-31']

def plot_march_sept_daily_evolution(df_kz, df_glob, title):
    # Process Adventdalen data
    df_kz['hour_min'] = df_kz.index.strftime('%H:%M')
    df_kz['month'] = df_kz.index.month
    df_kz = df_kz[df_kz['month'].between(4, 9)]
    monthly_avg_kz = df_kz.groupby(['month', 'hour_min']).mean().reset_index()
    pivot_kz = monthly_avg_kz.pivot(index='hour_min', columns='month', values='GHI')

    # Process GLOB data
    df_glob['hour_min'] = df_glob.index.strftime('%H:%M')
    df_glob['month'] = df_glob.index.month
    df_glob = df_glob[df_glob['month'].between(4, 9)]
    monthly_avg_glob = df_glob.groupby(['month', 'hour_min']).mean().reset_index()
    pivot_glob = monthly_avg_glob.pivot(index='hour_min', columns='month', values='GHI')

     # Define colors for each month
    colors = plt.cm.tab10(np.linspace(0, 1, 7))  # 7 colors for March-Sept

    # Plot
    plt.figure(figsize=(12, 6))
    for i, month in enumerate(range(4,10)):
        # Plot Adventdalen: solid line with circle markers every hour
        hourly_kz = pivot_kz.loc[[h for h in pivot_kz.index if h.endswith(':00')], month]
        plt.plot(pivot_kz.index, pivot_kz[month], '--',
                 color=colors[i], label=f'KZ ({month_abbr[month]})')
        plt.scatter(hourly_kz.index, hourly_kz, color=colors[i], marker='.', s=100, zorder=3)
        
        # Plot GLOB: dashed line with square markers every hour
        hourly_glob = pivot_glob.loc[[h for h in pivot_glob.index if h.endswith(':00')], month]
        plt.plot(pivot_glob.index, pivot_glob[month], '-',
                 color=colors[i], label=f'GLOB ({month_abbr[month]})')
        plt.scatter(hourly_glob.index, hourly_glob, color=colors[i], marker='.', s=100, zorder=3)
        
    # Add vertical line at 11:00
    plt.axvline(x=11, color='black', linewidth=2, label='11:00')

    # Customize plot
    # plt.title(f'{title}')
    plt.xlabel('Time of Day')
    plt.ylabel('GHI [W/m²]')
    full_hours = [f"{h:02d}:00" for h in range(24)]
    plt.xticks(full_hours, full_hours, rotation=45)
    plt.legend(title='Dataset & Month', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid()
    plt.tight_layout()
    plt.show()

# Call the function
plot_march_sept_daily_evolution(df_kz, df_glob, 'GHI Daily Evolution')