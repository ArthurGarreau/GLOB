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
from itertools import combinations
from datetime import datetime
from config_path import SCRIPT_PATH, DATA_PATH
import glob_functions_calculation as fct

method = 'linear'; year = 2025; f = 10 #min data frequency

############################## File Paths #####################################

bsrn_datafile = DATA_PATH.parent / "Irradiance" / "NYA" / "NYA_radiation_2025-all.tab"
input_file = DATA_PATH / f"GLOB_data_{f}min_{year}.nc"

output_file_path = SCRIPT_PATH / "Data" / "Beam_Diffuse_Estimations"
output_file_bestestim = output_file_path / \
    f"{year}_bestestimation_beam_diffuse_{f}min_{method}.csv"

###############################################################################


# Define the azimuth directions and angles of the pyranometers to use from GLOB
azimuth_directions = ['N', 'E', 'S', 'W']; angles = [45, 90, 135]
pyrano_var = ['GHI'] + [f"{azimuth}_{angle}" for azimuth in azimuth_directions for angle in angles]
# Generate all combinations of pyrano_var
combs = list(combinations(pyrano_var, 3))

# Load NYA data (true_estimation)
df_NYA = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')
df_NYA.index = df_NYA.index.tz_localize('UTC')
df_NYA = df_NYA.resample(f'{f}min').first()
true_estimations = df_NYA[['DIF', 'DIR']]

# Load GLOB data
ds_glob = xr.open_dataset(input_file)
lat_glob = float(ds_glob.latitude.values); lon_glob = float(ds_glob.longitude.values)

# List to store all results
all_results = []
# Create a daily date range for the specified month and year
if year == 2023: start_date, end_date = (f'{year}-03-14', f'{year}-10-13')
if year == 2024: start_date, end_date = (f'{year}-04-14', f'{year}-10-13')
if year == 2025: start_date, end_date = (f'{year}-03-16', f'{year}-03-31')
dates = pd.date_range(start=start_date, end=end_date, freq='D')
# % Beam and Diffuse Irradiance Calculation
for date in dates:
    date = date.strftime('%Y-%m-%d')

    df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()
    df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')

    # Loop over each minute in the day
    TIMESTAMPS = pd.date_range(
        start=f"{date} 00:00:00", end=f"{date} 23:50:00", freq=f'{f}min', tz='UTC')

    for timestamp in TIMESTAMPS:
        glob_value = df_glob_one_day.loc[timestamp]
        # Calculate the solar position
        solar_angles = fct.calculate_solar_angles(timestamp, lat_glob, lon_glob)

        # Calculate the estimation based on the least square error method
        true_estimation = np.array(true_estimations.loc[timestamp])
        new_result = fct.find_best_estimation(combs,
                                               glob_value,
                                               solar_angles,
                                               lat_glob, lon_glob,
                                               true_estimation,
                                               method=method)
        all_results.append([str(timestamp)] + new_result) # the pyranometers combination is returned

        print(str(timestamp))

# Prepare the storage of all the results
column_names = ['Timestamp', 'Diffuse', 'Beam', 'Albedo',
        'Diffuse_prime', 'Beam_prime', 'Pyrano_Combination']
df_all_results = pd.DataFrame(all_results, columns=column_names)
df_all_results.set_index('Timestamp', inplace=True)

most_used = fct.most_frequent_pyrano_combination(df_all_results['Pyrano_Combination'])
pyrano_var = str(most_used[0]) + " with frequency " + str(most_used[1]) + "%"

# Write the header and units to the file
# Get the current date
current_date = datetime.now().strftime("%Y-%m-%d")

header = \
f"\
# Date of production: {current_date}\n\
# Produced by: Arthur Garreau\n\
# Estimation of beam and diffuse irradiance with GLOB using the least square error from the Faiman et al. (1992) method.\n\
# In the linear method the albedo is measured, in the non linear the albedo is estimated.\n\
# Location: {lat_glob}N {lon_glob}E\n\
# The most used orientations are : {pyrano_var}\n\
#\n\
#\n\
#\n\
#[UTC]\t\t[W m-2]\t [W m-2]\t [/]\t [W m-2]\t [W m-2]t [/]\n"

output_file = output_file_bestestim

# Open the file and write the header and units
with open(output_file, 'w', encoding='utf-8') as file:
    file.write(header)

# Append the DataFrame to the file
df_all_results.to_csv(output_file, index=True, mode='a', sep='\t', na_rep='NaN', encoding='utf-8')

ds_glob.close()


