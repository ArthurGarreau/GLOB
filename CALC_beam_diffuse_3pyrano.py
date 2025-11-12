# -*- coding: utf-8 -*-
"""
Beam and Diffuse Irradiance Calculation for GLOB
================================================

This script is similar to CALC_beam_diffuse.py, but the estimations are made
with a selected set of 3 pyranometers among GLOB pyranometers. 

Key Features:
-------------
- Loads GLOB data from a NetCDF file.
- Estimates beam and diffuse irradiance for each timestamp with 3 pyranometers.
- Calculate the metric of beam and diffuse estimations for the year 2025.

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: May 30, 2025
"""

# %% Load Libraries
import pandas as pd
from datetime import datetime
from config_path import DATA_PATH
import glob_functions_calculation as fct

method = 'linear'; year = 2025; f = 10 #min data frequency

pyrano_vars = [['GHI', 'N_90', 'N_45'],
               ['GHI', 'E_45', 'N_45'],
               ['GHI', 'S_45', 'N_45'],
               ['GHI', 'E_45', 'W_45']]

for pyrano_var in pyrano_vars:
    ############################## File Paths #####################################
    glob_datafile = DATA_PATH / "GLOB_data"  / f"GLOB_data_{f}min_{year}.nc"
    
    output_file = DATA_PATH / "Estim_Beam_Diffuse"  \
     / f"{year}_estimation_beam_diffuse_{f}min_{method}_{(pyrano_var)}pyrano.csv"
    ###############################################################################
    
    # Load GLOB data
    ds_glob = fct.read_netcdf(glob_datafile)

    lat_glob = float(ds_glob.lat.values); lon_glob = float(ds_glob.lon.values)
    
    # List to store all results
    all_results = []
    # Create a daily date range for the specified month and year
    if year == 2023: start_date, end_date = (f'{year}-03-14', f'{year}-10-13')
    if year == 2024: start_date, end_date = (f'{year}-04-14', f'{year}-10-13')
    if year == 2025: start_date, end_date = (f'{year}-03-16', f'{year}-07-31')
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
            new_result = fct.estimate_diffuse_beam(pyrano_var,
                                                    glob_value,                                                            
                                                    lat_glob, lon_glob,
                                                    method=method)
            all_results.append([str(timestamp)] + new_result)
        print(date)
    
    # Prepare the storage of all the results
    column_names = ['Timestamp', 'Diffuse', 'Beam', 'Albedo',
                    'Diffuse_prime', 'Beam_prime']
    # Convert all results to a DataFrame at once
    df_all_results = pd.DataFrame(all_results, columns=column_names)
    df_all_results.set_index('Timestamp', inplace=True)
    
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
# The orientations used are : {pyrano_var}\n\
#\n\
#\n\
#\n\
#[UTC]\t\t[W m-2]\t [W m-2]\t [/]\t [W m-2]\t [W m-2]\n"
    
    # Create the folder if it doesn't exist already
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Open the file and write the header and units
    with open(output_file, 'w', encoding='utf-8') as file:
        file.write(header)
    # Append the DataFrame to the file
    df_all_results.to_csv(output_file, index=True, mode='a', sep='\t', na_rep='NaN', encoding='utf-8')
    
    ds_glob.close()

################################ END ##########################################