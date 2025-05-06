# -*- coding: utf-8 -*-
"""
Solar Irradiance Calculation Script
===================================

This script calculates the solar irradiance (W m-2) on multiple orientations using Eq. 2.1 from Faiman et al. (1992).
The calculations are performed for a specified range of tilt and azimuth angles over a given date range.

Key Features:
-------------
- Loads daily solar data from CSV files.
- Computes the incident and zenith angles for each combination of tilt and azimuth angles.
- Calculates the irradiance using the beam and diffuse components of solar radiation.
- Filters out negative irradiance values.
- Saves the results to CSV files with a detailed header including metadata.

Dependencies:
-------------
- numpy
- pandas
- matplotlib
- pathlib
- datetime
- sys
- os
- Custom module: glob_functions_Faiman

Author: Arthur Garreau
Date: April 24, 2025
"""
import numpy as np
import pandas as pd
import pvlib
from datetime import datetime, timedelta
import sys
import os
from config import SCRIPT_PATH, DATA_PATH

f = 5 #min data frequency

############################## File Paths #####################################

# Add script path for the importing the functions in glob_functions_Faiman.py
sys.path.append(str(SCRIPT_PATH))

input_file_path_LYR = SCRIPT_PATH / "B_and_D_Estimations_LYR"
output_file_path_LYR = SCRIPT_PATH / "MultiOrientations_Irradiance_Data_LYR"

###############################################################################

# % Longyearbyen 

import glob_functions_Faiman as fct

# Define the location
latitude = 78.200318
longitude = 15.840308

# Define the range of tilt and azimuth angles
tilt_angles = np.arange(0, 91, 10)
azimuth_angles_calc = np.arange(-180, 180, 15)
azimuth_angles_names = np.arange(0, 360, 15)

# Define the date range
start_date = datetime(2024, 10, 1)
end_date = datetime(2024, 10, 5)

# Loop through each day
current_date = start_date
while current_date <= end_date:
    # Generate the file path for the current date
    file_date_str = current_date.strftime("%Y-%m-%d")
    filename = input_file_path_LYR / f"best_estimations_{file_date_str}_ERBS.csv"

    # Check if the file exists
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        current_date += timedelta(days=1)
        continue

    # Load the data
    try:
        data = pd.read_csv(filename, parse_dates=['Timestamp'], index_col='Timestamp', sep='\t', header=5)
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
    
    output_file_path_LYR.mkdir(parents=True, exist_ok=True)
    # Save the results to a CSV file for the current day
    output_file =  output_file_path_LYR / f"irradiance_results_{file_date_str}.csv"
    
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
