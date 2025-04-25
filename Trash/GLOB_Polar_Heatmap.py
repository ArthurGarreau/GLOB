# -*- coding: utf-8 -*-
"""
Polar Heatmap Generation Script for Solar Irradiance Data
=========================================================

This script generates polar heatmaps to visualize the average solar irradiance for multiple orientations.
It processes both monthly and daily data, calculating the average irradiance for specified tilt and azimuth angles,
and saves the resulting heatmaps as images.

Key Features:
-------------
- Loads solar irradiance data from CSV files.
- Filters data based on specified time ranges.
- Calculates the average irradiance for each combination of tilt and azimuth angles.
- Generates polar heatmaps using Matplotlib.
- Saves the heatmaps as PNG files with detailed titles and labels.

Dependencies:
-------------
- matplotlib
- numpy
- pandas
- pathlib
- datetime
- os

Author: Arthur Garreau
Date: April 25, 2025
"""

# %% Load Libraries
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


############################## File Paths #####################################

data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB")
output_file_path = data_path / "B_and_D_Estimations_LYR" 
save_fig_path = Path(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\Fig")

###############################################################################

# Load GLOB Data
ds_glob = xr.open_dataset(data_path / "GLOB_data_5min_2023-24.nc")

# %% Plot 

# Assuming ds_glob is your xarray Dataset

# Define the azimuth orientations and their corresponding angles in degrees
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
azimuth_angles = np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315, 360])

# Define the inclination angles
inclination_angles = [0, 45, 90, 135]

# Define the time range for filtering
start_date = '2024-03-01'
end_date = '2024-06-01'
start_time = 0
end_time = 23
period = 'Spring'
year = start_date[0:4]

# Filter the dataset by date and time
filtered_ds = ds_glob.sel(Timestamp=slice(start_date, end_date))
filtered_ds = filtered_ds.sel(Timestamp=(filtered_ds['Timestamp'].dt.hour >= start_time) & (filtered_ds['Timestamp'].dt.hour <= end_time))

# Create a grid for the heatmap
theta, r = np.meshgrid(azimuth_angles, inclination_angles)
irradiance_grid = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)

# Populate the grid with mean irradiance values
for i, inclination in enumerate(inclination_angles):
    for j, azimuth in enumerate(azimuth_orientations):
        if inclination == 0 and azimuth == 'N':
            # Use GHI for azimuth = 0 and angle = 0
            var_name = 'GHI'
        else:
            var_name = f"{azimuth}_{inclination}"



        if var_name == 'GHI':
            mean_irradiance = filtered_ds[var_name].mean().item()
            irradiance_grid[i,:] = mean_irradiance * np.ones(9)
        else:
            if var_name in filtered_ds:
                mean_irradiance = filtered_ds[var_name].mean().item()
                irradiance_grid[i, j] = mean_irradiance

# Duplicate the first column to close the loop in the polar plot
# irradiance_grid = np.hstack([irradiance_grid, irradiance_grid[:, :1]])
# theta = np.hstack([theta, theta[:, :1]])  # Duplicate the first column in theta as well
# r = np.hstack([r, r[:, :1]])  # Duplicate the first column in r as well

# Create a polar heatmap with contour lines
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

# Set the zero direction to point up (North) and rotate the plot
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

# Plot the heatmap with contour lines using the 'jet' colormap
contour = ax.contourf(theta, r, irradiance_grid, cmap='jet', levels=np.linspace(0, 100*np.round(np.max((irradiance_grid+50)/100)), 10), alpha=1)
contour_lines = ax.contour(theta, r, irradiance_grid, colors='black', levels=contour.levels, linewidths=0.5)
# Add contour line labels
# ax.clabel(contour_lines, fontsize=8, fmt='%d')

plt.colorbar(contour, label='Irradiance (W/m²)')

# Set the labels and title with fontsize 12
ax.set_xticks(azimuth_angles)
ax.set_xticklabels(azimuth_orientations, fontsize=12)
ax.set_yticks(inclination_angles)
ax.set_yticklabels(inclination_angles, fontsize=12)
ax.spines['polar'].set_color('black')
ax.tick_params(axis='both', colors='black')
ax.set_title(f'GLOB Polar Heatmap - {period} {year} - {start_time}:00-{end_time}:00 UTC', fontsize=12)

# Save the plot
save_fig_path.mkdir(parents=True, exist_ok=True)
fig.savefig(save_fig_path / f"polar_heatmap_{period}_2023_{start_time}-{end_time}UTC.png", dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

# %% Test timestamp


import matplotlib.pyplot as plt
import xarray as xr

# Assuming ds_glob is your xarray Dataset
# ds_glob = xr.open_dataset('your_dataset.nc')  # Uncomment and adjust if loading from a file

# Define the date for filtering
date_to_plot = '2024-05-02'

# Filter the dataset for the specific date
ghi_data = ds_glob['GHI'].sel(Timestamp=date_to_plot)
s_45_data = ds_glob['S_45'].sel(Timestamp=date_to_plot)
s_90_data = ds_glob['S_90'].sel(Timestamp=date_to_plot)

# Plot the data
plt.figure(figsize=(12, 8))

# Plot GHI
ghi_data.plot(label='GHI')

# Plot S_45
s_45_data.plot(label='S_45')

# Plot S_90
s_90_data.plot(label='S_90')

# Add title and labels
plt.title('Irradiance Components on 02-05-2024')
plt.xlabel('Time')
plt.ylabel('Irradiance (W/m²)')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()


