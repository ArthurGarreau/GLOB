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


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import calendar
from datetime import datetime
import os
from config_path import SCRIPT_PATH, DATA_PATH

method = 'nonlinear'; year = 2025; f = 10 #min data frequency
pyrano_var = np.zeros(25)

############################## File Paths #####################################

data_gti_dir_path = SCRIPT_PATH / "Data" / "GTI_Estimations"
gti_estimation_datafile = SCRIPT_PATH / 'Data' / "GTI_Estimations" / \
    f"{year}_estimation_GTI_{f}min_{method}_{len(pyrano_var)}pyrano.csv"
    
output_monthly_dir = SCRIPT_PATH / "Figures_all" / "Polar_Heatmap" / "Monthly"
###############################################################################


# %% Polar Heat Map - Monthly and daily average


gti_data = pd.read_csv(gti_estimation_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)
    
# Define the date range
if year == 2023: start_date, end_date = (f'{year}-03-14', f'{year}-10-13') 
if year == 2024: start_date, end_date = (f'{year}-04-14', f'{year}-10-13') 
if year == 2025: start_date, end_date = (f'{year}-03-16', f'{year}-03-31') 

date_range = pd.date_range(start_date, end_date)
first_month = date_range[0].month; last_month = date_range[-1].month
monthly_date_ranges = [
    pd.date_range(
        start=f'2024-{month:02d}-01',
        end=f'2024-{month:02d}-{pd.Timestamp(f"2024-{month:02d}-01").days_in_month}',
        freq='D'
    )
    for month in range(first_month, last_month)
]
# Define the time range for filtering
start_hour = 0
end_hour = 23

# Define the azimuth orientations and their corresponding angles in degrees
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
azimuth_angles = np.array([0, 45, 90, 135, 180, 225, 270, 315, 360])

# Define the inclination angles
inclination_angles = np.arange(0, 100, 10)
# inclination_angles = np.arange(90, 180, 10)

# Initialize a grid to accumulate irradiance values for monthly averages
theta, r = np.meshgrid(np.deg2rad(azimuth_angles), inclination_angles)

# Loop through each month's date range
for monthly_date_range in monthly_date_ranges:
    # Initialize accumulators for the current month
    irradiance_avg = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)
    day_count = 0

    # Get the current month
    current_month = monthly_date_range[0].month

    # Loop through each day in the current month's date range
    for date in monthly_date_range:
        try:
            # Attempt to filter the data for the current date
            daily_gti = gti_data.loc[str(date.date())]
        
        
        except KeyError:
            # If the date is not found, print a message and continue to the next date
            print(f"No data available for {date.date()}")
            continue
        
        # Filter the data between start_hour and end_hour
        filtered_data = daily_gti[(daily_gti.index.hour >= start_hour) & (daily_gti.index.hour < end_hour)]

        # Create daily grid
        daily_grid = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)

        # Populate the grid
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                var_name = f"gti{azimuth}_{inclination}"
                if var_name in filtered_data.columns:
                    mean_irradiance = np.nanmean(filtered_data[var_name])
                    daily_grid[i, j] = mean_irradiance

        # Duplicate first column to last for continuity
        daily_grid[:, -1] = daily_grid[:, 0]

        # Add to accumulator for monthly averages
        irradiance_avg = np.nansum([irradiance_avg, daily_grid], axis=0)
        day_count += 1

    # Calculate the monthly average
    irradiance_avg = irradiance_avg / day_count

    # Create the monthly polar heatmap
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # Plot with improved color scaling
    vmax = np.nanmax(irradiance_avg)
    contour = ax.contourf(theta, r, irradiance_avg, cmap='gnuplot2', levels=np.arange(0, 360, 10), alpha=1)

    # Add colorbar and formatting
    cbar = plt.colorbar(contour, label='Irradiance (W/m²)')
    ax.set_xticks(np.radians(azimuth_angles))
    ax.set_xticklabels(azimuth_orientations, fontsize=12)
    ax.set_yticks(inclination_angles)
    ax.set_yticklabels(inclination_angles, fontsize=12)
    ax.set_rlabel_position(150)
    ax.spines['polar'].set_color('black')
    ax.tick_params(axis='both', colors='black')
    ax.text(np.radians(150), inclination_angles[-1]+10, "Tilt angle (°)")
    ax.set_title(f'Monthly Average GLOB Polar Heatmap ({calendar.month_name[current_month]} 2024)\n{start_hour}:00-{end_hour}:00 UTC', fontsize=12)
    
    # plt.show()
    # Save and show the monthly plot
    output_monthly_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_monthly_dir / f"monthly_avg_polar_heatmap_{current_month}-2024.png", dpi=300, bbox_inches='tight')
    plt.close()
    


# %% Polar Heat Map - Seasonal average with BSRN data

year = 2025

data_daily_dir_path = data_nya_daily_dir_path
save_fig_monthly_path = save_fig_nya_monthly_path
save_fig_daily_path = save_fig_nya_daily_path

    
# Define the date range
start_date = datetime(year, 3, 30)
end_date = datetime(year, 3, 30)
date_range = pd.date_range(start_date, end_date)

# Define the time range for filtering
start_hour = 0
end_hour = 23

# Define the azimuth orientations and their corresponding angles in degrees
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
azimuth_angles = np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315, 360])

# Define the inclination angles
inclination_angles = np.arange(0, 91, 10)

# Initialize a grid to accumulate irradiance values for monthly averages
irradiance_avg = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)
theta, r = np.meshgrid(azimuth_angles, inclination_angles)
day_count = np.zeros((10,9))+np.nan
current_month = start_date.month

# Loop through each day
for date in date_range:
    file_date_str = date.strftime("%Y-%m-%d")
    file_path = data_nyabsrn_daily_dir_path / f"irradiance_results_{file_date_str}.csv"

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    try:
        # Load the data
        data = pd.read_csv(file_path, header=5)
        data['Timestamp'] = pd.to_datetime(data['Timestamp'])

        # Filter the data between 10:00 and 14:00
        filtered_data = data[(data['Timestamp'].dt.hour >= start_hour) & (data['Timestamp'].dt.hour < end_hour)]

        # Create daily grid
        daily_grid = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)

        # Populate the grid
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                azimuth_deg = int(np.rad2deg(azimuth))
                var_name = f"R_{inclination}_{azimuth_deg}"
                if var_name in filtered_data.columns:
                    mean_irradiance = np.nanmean(filtered_data[var_name])
                    daily_grid[i, j] = mean_irradiance

        # Duplicate first column to last for continuity
        daily_grid[:, -1] = daily_grid[:, 0]

        # Add to accumulator for monthly averages
        irradiance_avg = np.nansum([irradiance_avg,daily_grid], axis=0) 
        day_count = np.nansum([day_count, np.ones((10,9))+daily_grid*0], axis=0)

        # Check if the month has changed
        if (date + pd.Timedelta(days=1)).month != current_month:
            irradiance_avg = irradiance_avg / day_count

            # Create the monthly polar heatmap
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)

            # Plot with improved color scaling
            vmax = np.nanmax(irradiance_avg)
            contour = ax.contourf(theta, r, irradiance_avg, cmap='jet',
                                 levels=np.linspace(3, 303, 101), alpha=1)

            # Add colorbar and formatting
            cbar = plt.colorbar(contour, label='Irradiance (W/m²)')
            ax.set_xticks(azimuth_angles)
            ax.set_xticklabels(azimuth_orientations, fontsize=12)
            ax.set_yticks(inclination_angles)
            ax.set_yticklabels(inclination_angles, fontsize=12)
            ax.spines['polar'].set_color('black')
            ax.tick_params(axis='both', colors='black')
            ax.text(270.4, 100, "Tilt angle (°)")
            ax.set_title(f'Monthly Average BSRN Polar Heatmap ({calendar.month_name[current_month]} 2024)\n{start_hour}:00-{end_hour}:00 UTC', fontsize=12)

            # Save and show the monthly plot
            save_fig_monthly_path.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_fig_monthly_path / f"monthly_avg_polar_heatmap_{current_month}-2024_BSRN.png", dpi=300, bbox_inches='tight')
            # plt.close()
            # Reset the accumulator for the new month
            irradiance_avg = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)
            day_count = np.zeros((10,9))+np.nan
            current_month = (date + pd.Timedelta(days=1)).month

        # Create the daily polar heatmap
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        # Plot with improved color scaling
        vmax = np.nanmax(daily_grid)
        contour = ax.contourf(theta, r, daily_grid, cmap='jet', 
                             levels=np.linspace(6, 606, 101), alpha=1)

        # Add colorbar and formatting
        cbar = plt.colorbar(contour, label='Irradiance (W/m²)')
        ax.set_xticks(azimuth_angles)  # Exclude the duplicated angle
        ax.set_xticklabels(azimuth_orientations, fontsize=12)
        ax.set_yticks(inclination_angles)
        ax.set_yticklabels(inclination_angles, fontsize=12)
        ax.spines['polar'].set_color('black')
        ax.tick_params(axis='both', colors='black')
        ax.text(270.4, 100, "Tilt angle (°)")
        ax.set_title(f'Daily Average BSRN Polar Heatmap ({file_date_str})\n{start_hour}:00-{end_hour}:00 UTC', fontsize=12)

        # Save and show the daily plot
        save_fig_daily_path.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_fig_daily_path / f"daily_avg_polar_heatmap_{file_date_str}_BSRN.png", dpi=300, bbox_inches='tight')
        plt.close()

    except Exception as e:
        print(f"Error processing {file_date_str}: {str(e)}")
        continue


# %% Polar heatmap difference with Ny-Ålesund 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import calendar
from datetime import datetime
import os
from config_path import SCRIPT_PATH

method = 'nonlinear'; year = 2025; f = 10 #min data frequency
pyrano_var = np.zeros(25)

############################## File Paths #####################################

data_gti_dir_path = SCRIPT_PATH / "Data" / "GTI_Estimations"

gti_glob_datafile = SCRIPT_PATH / 'Data' / "GTI_Estimations" / \
    f"{year}_estimation_GTI_{f}min_{method}_{len(pyrano_var)}pyrano.csv"
gti_nya_datafile = SCRIPT_PATH / 'Data' / "GTI_Estimations" / \
    f"{year}_estimation_GTI_{f}min_NYA.csv"
      
output_monthly_dir = SCRIPT_PATH / "Figures_all" / "Polar_Heatmap" / "Monthly"
###############################################################################


gti_glob_data = pd.read_csv(gti_glob_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)
gti_nya_data = pd.read_csv(gti_nya_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)

gti_glob_data=gti_glob_data.resample('1H').mean()
gti_nya_data=gti_nya_data.resample('1H').mean()

gti_difference = abs(gti_glob_data-gti_nya_data)/(gti_glob_data+gti_nya_data)*2*100

gti_difference = gti_difference.resample('1M').mean()

irradiance_avg_mat = gti_difference.values.reshape(24, 37).transpose()
selected_rows = irradiance_avg_mat[0:36:2, :]; 
irradiance_avg_mat = selected_rows[:, 0:24:3]


# Define the azimuth orientations and their corresponding angles in degrees
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
azimuth_angles = np.array([0, 45, 90, 135, 180, 225, 270, 315, 360])

# Define the inclination angles
inclination_angles = np.arange(0, 180, 10)

# Initialize a grid to accumulate irradiance values for monthly averages
theta, r = np.meshgrid(np.deg2rad(azimuth_angles), inclination_angles)

# Duplicate first column to last for continuity
irradiance_avg_mat = np.hstack((irradiance_avg_mat, irradiance_avg_mat[:,[0]]))

# Create the monthly polar heatmap
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
# Plot with improved color scaling
vmax = np.nanmax(irradiance_avg_mat)
contour = ax.contourf(theta, r, irradiance_avg_mat, cmap='Blues', levels=np.arange(0, 80, 4), alpha=0.8)

# Add colorbar and formatting
cbar = plt.colorbar(contour, label='%')
ax.set_xticks(np.radians(azimuth_angles))
ax.set_xticklabels(azimuth_orientations, fontsize=12)
ax.set_yticks(inclination_angles[::2])
ax.set_yticklabels(inclination_angles[::2], fontsize=11)
ax.set_rlabel_position(150)
# ax.spines['polar'].set_color('green')
ax.tick_params(axis='both', colors='black')
ax.text(np.radians(150), inclination_angles[-1]+10, "Tilt angle (°)", fontsize=11)
ax.set_title('RPD :' + r' $GTI_{GLOB}$' + f' method {method}, {len(pyrano_var)} pyranometers', loc='left', fontsize=12)
# plt.subplots_adjust(left=0.1, right=0.5, top=0.7, bottom=0.1)

plt.show()
# Save and show the monthly plot
output_monthly_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(output_monthly_dir / f"monthly_avg_polar_heatmap_{method}_{year}_MBE_NYA_GLOB.png", dpi=300, bbox_inches='tight')
plt.close()




