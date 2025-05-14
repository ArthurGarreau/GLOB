# -*- coding: utf-8 -*-
"""
Beam and Diffuse Irradiance Calculation for GLOB in Longyearbyen
================================================================

This script calculates beam and diffuse irradiance using the Faiman et al. (1992) method for GLOB data in Longyearbyen.
It processes GLOB data (GLOB_data_5min_2023-24.nc), estimates irradiance components, generates plots, and compares the results with reference data.

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
import sys
from itertools import combinations
from datetime import datetime
from config import SCRIPT_PATH, DATA_PATH

f = 10 #min data frequency

############################## File Paths #####################################

# Add script path for the importing the functions in glob_functions_Faiman.py
sys.path.append(str(SCRIPT_PATH))

input_file = DATA_PATH / f"GLOB_data_{f}min_2023-24.nc"
output_file_path = SCRIPT_PATH / "Data" / "B_and_D_Estimations_LYR" 

###############################################################################

# Load GLOB data
ds_glob = xr.open_dataset(input_file)
lat_glob = float(ds_glob.latitude.values)
lon_glob = float(ds_glob.longitude.values)


# % Beam and Diffuse Irradiance Calculation
import glob_functions_Faiman_comb as fct


year = 2024
# Create a daily date range for the specified month and year
start_date = f'{year}-05-03'
end_date = f'{year}-05-05'
# end_date = f'{year}-{month:02d}-{pd.Period(start_date).days_in_month}'
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Define the azimuth directions and angles
azimuth_directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
angles = [45,90]

# Variables of interest for irradiance calculation
variables =  [f"{azimuth}_{angle}" for azimuth in azimuth_directions for angle in angles]
# Generate all combinations of variables
combs = list(combinations(variables, 2))
for date in dates: 
    date=date.strftime('%Y-%m-%d')
    
    df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()
    df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')
        
    # Prepare to store the results
    results = []
    
    # Loop over each minute in the day
    TIMESTAMPS = pd.date_range(start=f"{date} 00:00:00", end=f"{date} 23:10:00", freq=f'{f}min', tz='UTC')
    for timestamp in TIMESTAMPS:
        timestamp = pd.DatetimeIndex([timestamp])  # Ensure timestamp is in the desired format
        glob_value = df_glob_one_day.loc[timestamp]
        
        # Calculate the solar position
        solar_angles = fct.calculate_solar_angles(timestamp, lat_glob, lon_glob)
 
        # Use the optimized function
        # D, I, D_prime, I_prime, D_prime_var, D_prime_var, comb_opt = fct.estimation_diffuse_beam(variables, glob_value, solar_angles, lat_glob, lon_glob)
        D, I, D_prime, I_prime, D_prime_var, D_prime_var, comb_opt = fct.find_best_combination(combs, glob_value, solar_angles, lat_glob, lon_glob)

        # Append the results for this timestamp
        results.append({
            'Timestamp': timestamp[0],
            'Beam': I,  # Beam value
            'Diffuse': D,  # Diffuse value
            'Beam_prime': I_prime,  # Beam value
            'Diffuse_prime': D_prime,  # Diffuse value
            'Best_Combination': str(comb_opt), # Pyrano combination used for best_estim
            'Solar zenith': round(solar_angles['zenith'].values[0], 2),
            'Albedo': round(glob_value['albedo'].values[0], 2)
        })

        print(f"Timestamp: {timestamp[0].strftime('%Y-%m-%d %H:%M:%S')} | I={I} and D={D} W/m2 | Combination:{comb_opt}") 
    
    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    
    # Write the header and units to the file
    # Get the current date
    current_date = datetime.now().strftime("%Y-%m-%d")
    # Define your name
    your_name = "Arthur Garreau"
    date_production = f"Date of production: {current_date}\n"
    author = f"Produced by: {your_name}\n"
    header = "Best estimation of beam and diffuse irradiance with GLOB using the Faiman et al. (1992) method.\n\
Location: Adventdalen (78.200318N 15.840308E).\n"
    units = "[UTC]\t\t[W m-2]\t [W m-2]\t  [W m-2]\t [W m-2]\t [/]\t [°]\t [/]\n"
    
    output_file_path.mkdir(parents=True, exist_ok=True)
    output_file = output_file_path / f"estimation_beam_diffuse_{date}_{f}min.csv"
    # output_file = output_file_path / f"estimation_beam_diffuse_{date}_bestcombination_{f}min.csv"

    # Open the file and write the header and units
    with open(output_file, 'w', encoding='utf-8') as file:
        file.write(date_production)
        file.write(author)
        file.write(header)
        file.write(units)
    
    # Append the DataFrame to the file
    results_df.to_csv(output_file, index=False, mode='a', sep='\t', na_rep='NaN', encoding='utf-8')


