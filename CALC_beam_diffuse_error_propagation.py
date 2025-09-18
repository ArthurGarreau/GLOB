# -*- coding: utf-8 -*-
"""
Beam and Diffuse Error Propagation with a Monte Carlo Approach
==============================================================

This script uses a Monte Carlo approach to understand how the measurement error 
propagates to the estimations error of beam and diffuse through the model.
It makes estimations of beam and diffuse by inputing a random noise  
to the input of the model, i.e. GLOB pyranometers data and albedo data.

Key Features:
-------------
- Loads GLOB data from a NetCDF file.
- Uses a Monte Carlo approach to calculate 1000 times beam and diffuse with 
  the model and output the error due to the input noise.
  
Dependencies:
-------------
- xarray
- pandas
- numpy
- scipy.stats
- datetime
- Custom module: glob_functions_calculation

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: July 29, 2024
"""

# %% Load Libraries
import xarray as xr
import pandas as pd
import numpy as np
import scipy.stats as stats
from datetime import datetime
from config_path import B_D_DATA_PATH, GLOB_DATA_PATH
import glob_functions_calculation as fct

### Modified function for performing Monte Carlo calculation


method = 'linear'; year = 2024; f = 10 #min data frequency
n_simulations = 500
error = 0.1

pyrano_vars = [['GHI', 'N_45', 'E_45', 'S_45', 'W_45'],
               ['GHI', 'N_45', 'N_90', 'N_135', 'E_45', 'E_90', 'E_135',
                'S_45', 'S_90', 'S_135', 'W_45', 'W_90', 'W_135']]

# pyrano_vars = [['GHI', 'N_45', 'N_90', 'N_135', 'E_45', 'E_90', 'E_135',
#      'S_45', 'S_90', 'S_135', 'W_45', 'W_90', 'W_135']]
############################## File Paths #####################################

glob_datafile = GLOB_DATA_PATH / f"GLOB_data_{f}min_{year}.nc"

output_file_path = B_D_DATA_PATH

###############################################################################

for pyrano_var in pyrano_vars:

    output_file = output_file_path / f"{year}_MCestimation_beam_diffuse_{method}_{len(pyrano_var)}pyrano_noiseAlbedo_error_{error}.csv"
    
    # Load GLOB data
    ds_glob = xr.open_dataset(glob_datafile)
    lat_glob = float(ds_glob.latitude.values); lon_glob = float(ds_glob.longitude.values)
    
    # List to store all results
    all_results = []
    # Create a daily date range for the specified month and year
    start_date, end_date = (f'{year}-04-19', f'{year}-08-14')
    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    # % Beam and Diffuse Irradiance Calculation
    for date in dates:
        date = date.strftime('%Y-%m-%d')
    
        df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()
        df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')
    
        # Loop over each minute in the day
        TIMESTAMPS = pd.date_range(
            start=f"{date} 6:00:00", end=f"{date} 20:00:00", freq=f'1h', tz='UTC')
    
        for timestamp in TIMESTAMPS:
            glob_value = df_glob_one_day.loc[timestamp]
            # Calculate the solar position
            solar_angles = fct.calculate_solar_angles(timestamp, lat_glob, lon_glob)
            
            # Calculate the estimation based on the least square error method
            D_MC, B_MC = fct.estimation_diffuse_beam_MonteCarlo(pyrano_var,
                                                            glob_value,
                                                            solar_angles,
                                                            lat_glob, lon_glob,
                                                            n_simulations, error)
            
            error_D_MC, error_B_MC = np.nanstd(D_MC)/np.nanmean(D_MC), np.nanstd(B_MC)/np.nanmean(B_MC)
            stat, p_D_MC = stats.shapiro(D_MC)
            stat, p_B_MC = stats.shapiro(B_MC)
            
            all_results.append([str(timestamp)] + [error_D_MC, error_B_MC, p_D_MC, p_B_MC])
        print(date)
    
    # Prepare the storage of all the results
    column_names = ['Timestamp', 'Error Diffuse', 'Error Beam', "pvalue Diffuse", "pvalue Beam"]
    # Convert all results to a DataFrame at once
    df_all_results = pd.DataFrame(all_results, columns=column_names).round(3)
    df_all_results.set_index('Timestamp', inplace=True)
    
    avg_error_D_MC = df_all_results['Error Diffuse'].mean()*100
    avg_error_B_MC = df_all_results['Error Beam'].mean()*100
    percent_pvalue_D_MC = np.sum(df_all_results["pvalue Diffuse"] < 0.05) / len(df_all_results["pvalue Diffuse"]) *100
    percent_pvalue_B_MC = np.sum(df_all_results["pvalue Beam"] < 0.05) / len(df_all_results["pvalue Beam"]) *100
    
    # Write the header and units to the file
    # Get the current date
    current_date = datetime.now().strftime("%Y-%m-%d")
    
    header = \
    f"\
# Date of production: {current_date}\n\
# Produced by: Arthur Garreau\n\
# Error propagation simulation with a Monte Carlo approach of {n_simulations} simulations,\n\
# for the model of beam and diffuse estimations, with an input error of {round(error*100)}%. \n\
# The p-value correspond of a normality distribution test on the model output after each Monte Carlo simulatins.\n\
# Location: {lat_glob}N {lon_glob}E\n\
# The orientations used are : {pyrano_var}\n\
# GENERAL RESULTS:\n\
# Error model Diffuse = {avg_error_D_MC.round(2)} % \n\
# Error model Beam = {avg_error_B_MC.round(2)} % \n\
# % of p-value<0.05 Diffuse = {percent_pvalue_D_MC.round(2)} % \n\
# % of p-value<0.05 Beam = {percent_pvalue_B_MC.round(2)} % \n\
#[UTC]\t\t[%]\t [%]\t [%]\t [%]\n"
        
    # Create the folder if it doesn't exist already
    output_file.mkdir(parents=True, exist_ok=True)
    
    # Open the file and write the header and units
    with open(output_file, 'w', encoding='utf-8') as file:
        file.write(header)
    
    # Append the DataFrame to the file
    df_all_results.to_csv(output_file, index=True, mode='a', sep='\t', na_rep='NaN', encoding='utf-8')

    ds_glob.close()
