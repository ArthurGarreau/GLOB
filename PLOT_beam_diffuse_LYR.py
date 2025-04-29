# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 14:08:11 2025

@author: arthurg
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from datetime import datetime, timedelta

############################## File Paths #####################################

data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data")
kz_datafile_path = data_path / "Irradiance_ncdf" / "Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"
glob_estim_path = data_path / "GLOB" / "B_and_D_Estimations_LYR"
output_plot_path = Path(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\Fig\Estim_Faiman\lyr")

###############################################################################
ds_glob = xr.open_dataset(data_path / "GLOB" / "GLOB_data_5min_2023-24.nc")

# % Plot Beam, Diffuse, and Global irradiances estimated by GLOB
# Define criteria and date range
Criteria = "ERBS"
dates = pd.date_range(start="2024-05-02", end="2024-05-02", freq="D", tz="UTC")

for date in dates:
        
    date=date.strftime('%Y-%m-%d')
    results_df = pd.read_csv(glob_estim_path / f'best_estimations_{date}_{Criteria}_bis.csv', sep='\t', header=5)  # Adjust path as needed
    results_df['Timestamp'] = pd.to_datetime(results_df['Timestamp'])
    timestamps = results_df['Timestamp']
    
    GHI = ds_glob['GHI'].sel(Timestamp=date).to_pandas()
    zenith = np.deg2rad(results_df['Solar zenith'])
    
    plt.figure(figsize=(9, 6))
    #
    # Plot LYR Beam
    plt.plot(results_df['Timestamp'], results_df['Error'], "g:", label='Error')
    # Plot LYR Beam
    plt.plot(results_df['Timestamp'], results_df['Beam'], "r", label='Beam from GLOB')
    # Plot LYR Diffuse
    plt.plot(results_df['Timestamp'], results_df['Diffuse'], "r--", label='Diffuse from GLOB')
    # Plot LYR Global
    plt.plot(GHI, "k" , label='Global measured with Kipp&Zonen')
    plt.plot(results_df['Timestamp'], np.cos(zenith)*results_df['Beam']+results_df['Diffuse'], "b",
             label='Global = cos(z).Beam + Diffuse\nfrom GLOB')
  
    # Customize the plot
    plt.xlabel(f'Time (UTC)\n{date}', fontsize=14); plt.ylabel('Irradiance [W m²]', fontsize=14)
    plt.legend(fontsize=9)
    plt.grid(True, linestyle=':')
    # plt.tight_layout()

    # Set hourly ticks
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=3))  # Tick every hour
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))  # Format as HH:MM
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    # Show the plot
    # plt.savefig(output_plot_path / f"Estim_Faiman_{Criteria}_{date}")
    # plt.close()
    plt.show()
    
# %% Plot Global Vs. reconstructed Global

# Load K&Z data
ds_kz = xr.open_dataset(data_path / r"Irradiance_ncdf\Adventdalen_global_horizontal_irradiances_LW_SW_all.nc")


import numpy as np
import matplotlib.dates as mdates
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

dates = pd.date_range(start="2024-04-14", end="2024-07-30", freq="D", tz="UTC")
df_LYR_ghi = ds_kz["SWdown"].to_dataframe()
df_LYR_ghi = df_LYR_ghi.tz_localize('UTC')
df_LYR_ghi = df_LYR_ghi.sort_index()
# Lists to store data for plotting
all_timestamps = []
all_global_kz = []
all_global_construct = []

# Accumulate data for each date
for date in dates:
    date_str = date.strftime('%Y-%m-%d')
    results_df = pd.read_csv(data_path / "GLOB" / "Best_estimations_Faiman" / f'best_estimations_{date_str}_ERBS.csv')

    # Convert the Timestamp column to datetime
    results_df['Timestamp'] = pd.to_datetime(results_df['Timestamp'])

    # Select the relevant data from df_LYR_ghi
    timestamps = results_df['Timestamp']
    df_LYR_ghi_sel = df_LYR_ghi.loc[timestamps[0]:timestamps.iloc[-1]]
    df_LYR_ghi_sel = df_LYR_ghi_sel.rename(columns={"SWdown": "Global K&Z"})

    # Calculate the irradiance
    zenith_rad = np.deg2rad(results_df['Solar zenith (°)'])
    results_df["Constructed Global"] = np.cos(zenith_rad) * results_df['Beam'] + results_df['Diffuse']

    results_df = pd.merge(df_LYR_ghi_sel, results_df, how='inner', left_on='time', right_on='Timestamp')
    # Append data to lists
    all_global_kz.extend(results_df['Global K&Z'])
    all_global_construct.extend(results_df['Constructed Global'])


# Create a mask for non-NaN values
mask = ~np.isnan(all_global_kz) & ~np.isnan(all_global_construct)
all_global_kz = np.array(all_global_kz); all_global_construct = np.array(all_global_construct)
# Filter the data
all_global_kz = all_global_kz[mask]
all_global_construct = all_global_construct[mask]
# %
# Plot all data
plt.figure(figsize=(12, 8))
# Set font size for all elements
plt.rcParams.update({'font.size': 12})
# Create scatter plot
plt.scatter(all_global_kz, all_global_construct, marker='.', color='k', alpha=0.5)
# Plot the line x = y
plt.plot([min(all_global_kz), max(all_global_kz)], [min(all_global_kz), max(all_global_kz)], color='red', linestyle='--', 
         label='x = y')
# Use Ridge Regression for better numerical stability
model = Ridge(alpha=1e-3)  # alpha is the regularization strength
model.fit(all_global_kz.reshape(-1, 1), all_global_construct)
regression_line = model.predict(all_global_kz.reshape(-1, 1))
coefficients = model.coef_[0]
intercept = model.intercept_
r2 = r2_score(all_global_kz, regression_line)

# Plot the regression line
plt.plot(all_global_kz, regression_line, color='red', 
         label=f'Regression\n y = {coefficients:.2f}.x + {intercept:.2f} | $R^2$: {r2:.2f}')
# Show plot
plt.show()
plt.xlabel('Constructed Global from GLOB (W/m²)')
plt.ylabel('Global from Kipp&Zonen (W/m²)')
plt.title('Data from 2024')
plt.legend()
plt.grid(True, linestyle=':')
plt.show()
save_path = "C:/Users/arthurg/OneDrive - Universitetssenteret på Svalbard AS/Documents/UNIS_PhD/PAPER_2/PAPER_2_Data_Analysis/Fig/"
# plt.savefig(save_path + f"Regression_KZ_GLOB")
# plt.close()