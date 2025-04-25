# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 10:57:24 2025

@author: arthurg
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime, timedelta
import sys
import os

# Add custom script path
sys.path.append(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\GLOB_scripts")
import glob_functions_Faiman as fct
data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data")
# %%
# Define the location
latitude = 78.200318
longitude = 15.840308

# Define the range of tilt and azimuth angles
tilt_angles = np.arange(0, 91, 10)
azimuth_angles_calc = np.arange(-180, 180, 15)
azimuth_angles_names = np.arange(0, 360, 15)

# Define the date range
start_date = datetime(2024, 10, 1)
end_date = datetime(2024, 10, 31)

# Loop through each day
current_date = start_date
while current_date <= end_date:
    # Generate the file path for the current date
    file_date_str = current_date.strftime("%Y-%m-%d")
    file_path = fr"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB\B_and_D_Estimations_LYR\best_estimations_{file_date_str}_ERBS.csv"

    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        current_date += timedelta(days=1)
        continue

    # Load the data
    try:
        data = pd.read_csv(file_path, parse_dates=['Timestamp'], index_col='Timestamp', sep='\t', header=5)
    except Exception as e:
        print(f"Error loading data for {file_date_str}: {e}")
        current_date += timedelta(days=1)
        continue

    # Extract necessary columns
    beam_prime = data['Beam_prime'].values
    diffuse_prime = data['Diffuse_prime'].values
    solar_zenith = data['Solar zenith'].values
    albedo = data['Albedo'].values

    # Initialize a dictionary to store results for the current day
    results_dict = {}
    
    # Calculate irradiance for each combination of tilt and azimuth angles
    for tilt in tilt_angles:
        for i in range(len(azimuth_angles_calc)):
            theta_i, theta_z = fct.incident_and_zenith_angle(data.index, tilt, azimuth_angles_calc[i], latitude, longitude)
            b = np.cos(np.radians(theta_i)) + albedo * np.cos(np.radians(solar_zenith)) * (1 - np.cos(np.radians(theta_i))) / 2
            d = (1 + np.cos(np.radians(solar_zenith))) / 2 + albedo * (1 - np.cos(np.radians(solar_zenith))) / 2
            R = b * beam_prime + d * diffuse_prime
            results_dict[f'R_{tilt}_{azimuth_angles_names[i]}'] = R
    
    # Create a DataFrame from the dictionary
    results = pd.DataFrame(results_dict, index=data.index)

    # Filter out negative values
    results = round(results)
    results = results.where(results >= 0, np.nan)
    results['Timestamp'] = results.index
    # Reorder the columns to move Timestamp to the first column
    columns = ['Timestamp'] + [col for col in results.columns if col != 'Timestamp']
    results = results[columns]
    
    # Save the results to a CSV file for the current day
    output_file = data_path / "GLOB" / "MultiOrientations_Irradiance_Data_LYR" / f"irradiance_results_{file_date_str}.csv"
    
    # Write the header and units to the file
    # Get the current date
    today = datetime.now().strftime("%Y-%m-%d")
    # Define your name
    your_name = "Arthur Garreau"
    date_production = f"Date of production: {today}\n"
    author = f"Produced by: {your_name}\n"
    header = "Estimation of solar irradiance for multiple orientations (refer to Eq. 2.1 in Faiman et al. (1992)). Each column represents an \
estimation for a specific tilt and azimuth, denoted as R_{tilt}_{azimuth}.\n\
Location: Adventdalen (78.200318N 15.840308E).\n"
    units = "[UTC]\t[W m-2]\n"
    
    # Open the file and write the header and units
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(date_production)
        f.write(author)
        f.write(header)
        f.write(units)
    
    try:
        results.to_csv(output_file, index=False, mode='a', sep=',', encoding='utf-8')
        print(f"Results for {file_date_str} saved to {output_file}")
    except Exception as e:
        print(f"Error saving results for {file_date_str}: {e}")

    # Move to the next day
    current_date += timedelta(days=1)
