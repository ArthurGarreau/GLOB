# -*- coding: utf-8 -*-
"""
GLOB Data Visualization Script
==============================

This script processes GLOB data and generates visualizations, including polar heatmaps and time series plots,
to analyze irradiance components for different orientations and time periods.

Key Features:
-------------
- Loads GLOB data from a NetCDF file.
- Filters data based on specified date and time ranges.
- Generates polar heatmaps to visualize mean irradiance values for various azimuth and inclination angles.
- Plots time series data for specific dates to compare different irradiance components.

Dependencies:
-------------
- xarray
- pandas
- numpy
- matplotlib
- pathlib
- sys
- itertools

Author: Arthur Garreau
Date: March 28, 2025
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import calendar
from pathlib import Path
from datetime import datetime, timedelta
import os


# %% Polar Heat Map - Monthly average

### Define the base directory and date range
base_dir = r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB\MultiOrientations_Irradiance_Data_LYR"
for month in range(4,10):
    year = 2024
    # month = 7  # For example, October
    # Create a daily date range for the specified month and year
    start_date = f'{year}-{month:02d}-01'
    end_date = f'{year}-{month:02d}-{pd.Period(start_date).days_in_month}'
    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    
    start_hour = 0; end_hour = 23
    date_range = pd.date_range(start=start_date, end=end_date)
    
    # Define the azimuth orientations and their corresponding angles in degrees
    azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
    azimuth_angles = np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315, 360])
    
    # Define the inclination angles
    inclination_angles = np.arange(0, 91, 10)
    
    ### Initialize a grid to accumulate irradiance values
    irradiance_sum = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)
    theta, r = np.meshgrid(azimuth_angles, inclination_angles)
    day_count = 0
    
    for date in date_range:
        file_date_str = date.strftime("%Y-%m-%d")
        file_path = os.path.join(base_dir, f"irradiance_results_{file_date_str}.csv")
        
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue
        
        try:
            # Load the data
            data = pd.read_csv(file_path, header=5)
            data['Timestamp'] = pd.to_datetime(data['Timestamp'])
            
            # Filter the data between 10:00 and 14:00
            filtered_data = data[(data['Timestamp'].dt.hour >= start_hour) & (data['Timestamp'].dt.hour <= end_hour)]
            
            # Create daily grid
            daily_grid = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)
            
            # Populate the grid
            for i, inclination in enumerate(inclination_angles):
                for j, azimuth in enumerate(azimuth_orientations[:-1]):  # Exclude the duplicate N
                    var_name = f"R_{inclination}_{azimuth_angles[j] * 180 / np.pi:.0f}"
                    if var_name in filtered_data.columns:
                        mean_irradiance = np.nanmean(filtered_data[var_name])
                        daily_grid[i, j] = mean_irradiance
            
            # Duplicate first column to last for continuity
            daily_grid[:, -1] = daily_grid[:, 0]
            
            # Add to accumulator
            irradiance_sum = np.nansum(np.stack((irradiance_sum, daily_grid)), axis=0)
            day_count += 1
            
        except Exception as e:
            print(f"Error processing {file_date_str}: {str(e)}")
            continue
    print(day_count)
    # Calculate monthly average
    if day_count > 0:
        irradiance_avg = irradiance_sum / day_count
    else:
        raise ValueError("No valid days found for processing")
    
    ##################### Create the polar heatmap ################################
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    # Plot with improved color scaling
    vmax = np.nanmax(irradiance_avg)
    contour = ax.contourf(theta, r, irradiance_avg, cmap='jet', 
                         levels=np.linspace(0, 300, 101), alpha=1)
    
    # Add colorbar and formatting
    cbar = plt.colorbar(contour, label='Irradiance (W/m²)')
    ax.set_xticks(azimuth_angles)
    ax.set_xticklabels(azimuth_orientations, fontsize=12)
    ax.set_yticks(inclination_angles)
    ax.set_yticklabels(inclination_angles, fontsize=12)
    ax.spines['polar'].set_color('black')
    ax.tick_params(axis='both', colors='black')
    ax.set_title(f'Monthly Average GLOB Polar Heatmap ({calendar.month_name[month]} 2024)\n{start_hour}:00-{end_hour}:00 UTC', fontsize=12)
    
    ### Save and show
    plt.show()
    
    save_path = Path(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\Fig\Polar_Heatmap_LYR\Monthly")
    save_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path / f"monthly_avg_polar_heatmap_{month}-2024.png", dpi=300, bbox_inches='tight')
    plt.close()
    




# %% Polar Heat Map - Daily average

# Define the base directory
base_dir = r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB\Multiorientations_irradiance_data_LYR"

# Define the date range for April
start_date = datetime(2024, 4, 14)
end_date = datetime(2024, 9, 30)
date_range = pd.date_range(start_date, end_date)

# Define the time range for filtering
start_hour = 0
end_hour = 23

# Define the azimuth orientations and their corresponding angles in degrees
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
azimuth_angles = np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315, 360])

# Define the inclination angles
inclination_angles = np.arange(0, 91, 10)

# Loop through each day in April
for date in date_range:
    file_date_str = date.strftime("%Y-%m-%d")
    file_path = os.path.join(base_dir, f"irradiance_results_{file_date_str}.csv")

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    try:
        # Load the data
        data = pd.read_csv(file_path, header=5)  # Skip the header rows

        # Convert the Timestamp column to datetime
        data['Timestamp'] = pd.to_datetime(data['Timestamp'])

        # Filter the data between 10:00 and 14:00
        filtered_data = data[(data['Timestamp'].dt.hour >= start_hour) & (data['Timestamp'].dt.hour <= end_hour)]

        # Create a grid for the heatmap
        theta, r = np.meshgrid(azimuth_angles, inclination_angles)
        irradiance_grid = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)

        # Populate the grid with mean irradiance values using nanmean
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                azimuth_deg = int(np.rad2deg(azimuth))
                var_name = f"R_{inclination}_{azimuth_deg}"
                if var_name in filtered_data.columns:
                    mean_irradiance = np.nanmean(filtered_data[var_name])
                    irradiance_grid[i, j] = mean_irradiance

        # Duplicate the first column to the last column for continuity
        irradiance_grid[:, -1] = irradiance_grid[:, 0]

        # Create the polar heatmap
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        # Plot with improved color scaling
        vmax = np.nanmax(irradiance_grid)
        contour = ax.contourf(theta, r, irradiance_grid, cmap='jet',
                             levels=np.linspace(0, 600, 101), alpha=1)

        # Add colorbar and formatting
        cbar = plt.colorbar(contour, label='Irradiance (W/m²)')
        ax.set_xticks(azimuth_angles)  # Exclude the duplicated angle
        ax.set_xticklabels(azimuth_orientations, fontsize=12)
        ax.set_yticks(inclination_angles)
        ax.set_yticklabels(inclination_angles, fontsize=12)
        ax.spines['polar'].set_color('black')
        ax.tick_params(axis='both', colors='black')
        ax.set_title(f'Daily Average GLOB Polar Heatmap ({file_date_str})\n{start_hour}:00-{end_hour}:00 UTC', fontsize=12)
       
        # Show the plot
        plt.show()
        
        # Save the plot
        save_path = Path(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\Fig\Polar_Heatmap_LYR\Daily")
        save_path.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path / f"daily_avg_polar_heatmap_{file_date_str}.png", dpi=300, bbox_inches='tight')
        plt.close()

    except Exception as e:
        print(f"Error processing {file_date_str}: {str(e)}")
        continue




