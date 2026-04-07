# -*- coding: utf-8 -*-
"""
Beam and Diffuse Irradiance Calculation for GLOB
================================================

This script is similar to CALC_beam_diffuse.py but it tests all the possible 
combinations of 3 pyranometers to estimate beam and diffuse and keep the 
best one each timestep in comparison to Ny-Ålesund BSRN measurements.
Key Features:
-------------
- Loads GLOB data from a NetCDF file.
- Find the best 3 pyranometers set up to estimate beam and diffuse irradiance
  each timestamp.

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 22, 2024
"""

# %% Load Libraries
import pandas as pd
import numpy as np
import xarray as xr
from itertools import combinations
from datetime import datetime
from config_path import DATA_PATH
import glob_functions_calculation as fct


method = 'linear'; year = 2025; f = 10 #min data frequency

pyrano_vars = ['GHI','N_45', 'N_90', 'NE_45', 'NE_90', 'E_45', 'E_90',
     'SE_45', 'SE_90', 'S_45', 'S_90', 'SW_45', 'SW_90',
     'W_45', 'W_90', 'NW_45', 'NW_90']

############################## File Paths #####################################
bsrn_datafile = DATA_PATH / "NYA_BSRN_data" / "NYA_radiation_2025-all.tab"
glob_datafile = DATA_PATH / "GLOB_data" / f"GLOB_data_{f}min_{year}.nc"

output_file = DATA_PATH / "NYA_BSRN_data"  / \
    f"{year}_bestestimation_beam_diffuse_{f}min_{method}_2.csv"
###############################################################################

# Generate all combinations of pyrano_var
combs = list(combinations(pyrano_vars, 3))

# Load NYA data (true_estimation)
df_NYA = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')
df_NYA.index = df_NYA.index.tz_localize('UTC')
df_NYA = df_NYA.resample(f'{f}min').first()
true_estimations = df_NYA[['DIF', 'DIR']]

# Load GLOB data
ds_glob = xr.open_dataset(glob_datafile, engine="h5netcdf")
lat_glob = float(ds_glob.lat.values); lon_glob = float(ds_glob.lon.values)

# List to store all results
all_results = []
# Create a daily date range for the specified month and year
start_date, end_date = (f'{year}-08-01', f'{year}-09-05')
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

        # Calculate the estimation based on the least square error method
        true_estimation = np.array(true_estimations.loc[timestamp])
        new_result = fct.find_best_estimation(combs,
                                               glob_value,                                               
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

# Create the folder if it doesn't exist already
output_file.parent.mkdir(parents=True, exist_ok=True)

# Open the file and write the header and units
with open(output_file, 'w', encoding='utf-8') as file:
    file.write(header)

# Append the DataFrame to the file
df_all_results.to_csv(output_file, index=True, mode='a', sep='\t', na_rep='NaN', encoding='utf-8')
ds_glob.close()

################################ END ##########################################
