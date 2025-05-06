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
from pathlib import Path
from datetime import datetime
import os

############################## File Paths #####################################

data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB")
data_lyr_daily_dir_path = data_path / "MultiOrientations_Irradiance_Data_LYR"
data_nya_daily_dir_path = data_path / "MultiOrientations_Irradiance_Data_NYA"
data_nyabsrn_daily_dir_path = data_path / "MultiOrientations_Irradiance_Data_NYA_bsrn"


paper2_path =  Path(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis")
save_fig_lyr_monthly_path = paper2_path / "GLOB_scripts" / "Figures_all" / "Polar_Heatmap" / "LYR" / "monthly"
save_fig_lyr_daily_path = paper2_path / "GLOB_scripts" / "Figures_all" / "Polar_Heatmap" / "LYR" / "daily"
save_fig_nya_monthly_path = paper2_path / "GLOB_scripts" / "Figures_all" / "Polar_Heatmap" / "NYA" / "monthly"
save_fig_nya_daily_path = paper2_path / "GLOB_scripts" / "Figures_all" / "Polar_Heatmap" / "NYA" / "daily"


###############################################################################


# %% Polar Heat Map - Monthly and daily average

year = 2025

if year == 2025:
    data_daily_dir_path = data_nya_daily_dir_path
    save_fig_monthly_path = save_fig_nya_monthly_path
    save_fig_daily_path = save_fig_nya_daily_path
else:
    data_daily_dir_path = data_lyr_daily_dir_path
    save_fig_monthly_path = save_fig_lyr_monthly_path
    save_fig_daily_path = save_fig_lyr_daily_path
    
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
    file_path = data_daily_dir_path / f"irradiance_results_{file_date_str}.csv"

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
            ax.set_title(f'Monthly Average GLOB Polar Heatmap ({calendar.month_name[current_month]} 2024)\n{start_hour}:00-{end_hour}:00 UTC', fontsize=12)

            # Save and show the monthly plot
            save_fig_monthly_path.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_fig_monthly_path / f"monthly_avg_polar_heatmap_{current_month}-2024.png", dpi=300, bbox_inches='tight')
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
        ax.set_title(f'Daily Average GLOB Polar Heatmap ({file_date_str})\n{start_hour}:00-{end_hour}:00 UTC', fontsize=12)

        # Save and show the daily plot
        save_fig_daily_path.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_fig_daily_path / f"daily_avg_polar_heatmap_{file_date_str}.png", dpi=300, bbox_inches='tight')
        plt.close()

    except Exception as e:
        print(f"Error processing {file_date_str}: {str(e)}")
        continue


# %% Polar Heat Map - Monthly and daily average with BSRN data

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

