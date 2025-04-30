# -*- coding: utf-8 -*-
"""
Beam and Diffuse Irradiance Calculation for GLOB in Ny-Ålesund
==============================================================

This script calculates beam and diffuse irradiance using the Faiman et al. (1992) method for GLOB data in Ny-Ålesund.
It processes GLOB data (GLOB_data_30sec_2025_NYA.nc), estimates irradiance components, generates plots, and compares the results with reference data.

Key Features:
-------------
- Loads GLOB data from a NetCDF file.
- Calculates solar angles and irradiance components using custom functions.
- Estimates beam and diffuse irradiance for each timestamp.
- Generates plots for beam, diffuse, and global irradiance.
- Compares reconstructed global irradiance with reference data using regression analysis.

Dependencies:
-------------
- xarray
- pandas
- numpy
- matplotlib
- pathlib
- sys
- itertools
- datetime
- sklearn
- Custom module: glob_functions_Faiman

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 22, 2024
"""

# %% Load Libraries
import xarray as xr
import pandas as pd
import numpy as np
from pathlib import Path
import sys
from itertools import combinations
from datetime import datetime

############################## File Paths #####################################
# Add script path for the importing the functions in glob_functions_Faiman.py
sys.path.append(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\GLOB_scripts")

data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data")
output_file_path = data_path / "GLOB" / "B_and_D_Estimations_NYA" 

###############################################################################

# Load GLOB data
ds_glob = xr.open_dataset(data_path / r"GLOB\GLOB_data_5min_2025.nc")
lat_glob = ds_glob.latitude.values
lon_glob = ds_glob.longitude.values

# % Beam and Diffuse Irradiance Calculation
import glob_functions_Faiman as fct

# Define criteria and date range
Criteria = "ERBS"

month = 4  # For example, October
year = 2025

# Create a daily date range for the specified month and year
start_date = f'{year}-03-16'
end_date = f'{year}-04-27'
# end_date = f'{year}-{month:02d}-{pd.Period(start_date).days_in_month}'
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Define the azimuth directions and angles
azimuth_directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
angles = [45, 90, 135]
# Variables of interest for irradiance calculation
variables = ['GHI'] + [f"{azimuth}_{angle}" for azimuth in azimuth_directions for angle in angles]


# Generate all combinations of variables
combs = list(combinations(variables, 2))
for date in dates: 
    date=date.strftime('%Y-%m-%d')
    
    df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()
    df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')

    # Generate all possible combinations
    combs = list(combinations(variables, 2))
        
    # Prepare to store the results
    results = []
    
    # Loop over each minute in the day
    TIMESTAMPS = pd.date_range(start=f"{date} 00:00:00", end=f"{date} 23:50:00", freq='10min', tz='UTC')
    for timestamp in TIMESTAMPS:
        timestamp = pd.DatetimeIndex([timestamp])  # Ensure timestamp is in the desired format
        glob_value = df_glob_one_day.loc[timestamp]
        
        # Calculate the solar position
        solar_position = fct.calculate_solar_angles(timestamp, lat_glob, lon_glob)
        # Extract the zenith angle
        zenith_angle = solar_position['zenith'].values[0]
 
        # Use the optimized function
        D, I, D_prime, I_prime, error, comb_opt = fct.find_best_combination(combs, glob_value, zenith_angle, lat_glob, lon_glob)
        
        comb_opt = str(comb_opt)
        if not np.isnan(D) and not np.isnan(I):
            D = round(D); I = round(I)
            D_prime = round(D_prime); I_prime = round(I_prime)
            error = round(error)
        else: 
            D = str(D); I = str(I)
            D_prime = str(D_prime); I_prime = str(I_prime)
            error = str(np.nan)
        zenith_angle = round(zenith_angle, 2)
        
        # Append the results for this timestamp
        results.append({
            'Timestamp': timestamp[0],
            'Beam': I,  # Beam value
            'Diffuse': D,  # Diffuse value
            'Beam_prime': I_prime,  # Beam value
            'Diffuse_prime': D_prime,  # Diffuse value
            'Error': error, # error = ((I_ref - I)**2 + (D_ref - D)**2)**0.5
            'Best_Combination': comb_opt, # Pyrano combination used for best_estim
            'Solar zenith': zenith_angle,
            'Albedo': round(glob_value['albedo'].values[0], 2)
        })

        print(f"Timestamp: {timestamp[0].strftime('%Y-%m-%d %H:%M:%S')} | I={I} and D={D} W/m2 | Combination:{comb_opt}")

    
    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    output_file = output_file_path / f"best_estimations_{date}_{Criteria}.csv"
    
    # Write the header and units to the file
    # Get the current date
    current_date = datetime.now().strftime("%Y-%m-%d")
    # Define your name
    your_name = "Arthur Garreau"
    date_production = f"Date of production: {current_date}\n"
    author = f"Produced by: {your_name}\n"
    header = "Best estimation of beam and diffuse irradiance with GLOB using the Faiman et al. (1992) method.\n\
Location: Ny-Alesund (78.92240N 11.92174E).\n"
    units = "[UTC]\t\t[W m-2]\t [W m-2]\t  [W m-2]\t [W m-2]\t [/]\t [/]\t [°]\t [/]\n"
    
    # Open the file and write the header and units
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(date_production)
        f.write(author)
        f.write(header)
        f.write(units)
    
    # Append the DataFrame to the file
    results_df.to_csv(output_file, index=False, mode='a', sep='\t', encoding='utf-8')

ds_glob.close()
