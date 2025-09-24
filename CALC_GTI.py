# -*- coding: utf-8 -*-
"""
Global Tilted Irradiance Calculation Script
===========================================

This script calculates the solar irradiance (W m-2) on multiple orientations 
using Eq. 2.1 from Faiman et al. (1992).
The calculations are performed for a specified range of tilt and azimuth angles
 over a given date range.

Key Features:
-------------
- Loads daily solar data from CSV files.
- Computes the incident and zenith angles for each combination of tilt and azimuth angles.
- Calculates the irradiance using the beam and diffuse components of solar radiation.
- Saves the results to CSV files with a detailed header including metadata.

Dependencies:
-------------
- numpy
- pandas
- datetime
- re
- pvlib
- Custom module: glob_functions_calculation

Author: Arthur Garreau
Date: April 24, 2025
"""

import numpy as np
import pandas as pd
import pvlib
from datetime import datetime
import re
from config_path import B_D_DATA_PATH, GTI_DATA_PATH
import glob_functions_calculation as fct

method = 'linear'; year = 2025; f = 10 #min data frequency

pyrano_vars = [
    ['GHI', 'N_45', 'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45', 'E_90', 'E_135',
     'SE_45', 'SE_90', 'SE_135', 'S_45', 'S_90', 'S_135', 'SW_45', 'SW_90', 'SW_135',
     'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135'],
    ['GHI', 'N_45', 'N_90', 'N_135', 'E_45', 'E_90', 'E_135',
     'S_45', 'S_90', 'S_135', 'W_45', 'W_90', 'W_135'],
    ['GHI', 'N_45', 'N_90', 'E_45', 'E_90', 'S_45', 'S_90', 'W_45', 'W_90'],
    ['GHI', 'N_45', 'E_45', 'S_45', 'W_45']
    ]
# pyrano_vars = [['GHI', 'N_90', 'N_45']]
# pyrano_vars = [['GHI', 'E_45', 'N_45']]


### Define the date range with timezone-aware datetime objects
if year == 2023: start_date, end_date = (f'{year}-03-14', f'{year}-10-13') 
if year == 2024: start_date, end_date = (f'{year}-04-14', f'{year}-10-13') 
if year == 2025: start_date, end_date = (f'{year}-03-16', f'{year}-07-31') 

### Define the range of tilt and azimuth angles
tilt_angles_calc = np.arange(0, 181, 5)
azimuth_angles_calc = np.arange(-180, 180, 15)
azimuth_angles = np.arange(0, 360, 15)

for pyrano_var in pyrano_vars:
    if len(pyrano_var)==3: pyrano_var_name = pyrano_var
    else: pyrano_var_name = len(pyrano_var)
        
    ############################## File Paths #####################################
    
    input_filename = B_D_DATA_PATH / \
        f"{year}_estimation_beam_diffuse_{f}min_{method}_{pyrano_var_name}pyrano.csv"
    output_file =  GTI_DATA_PATH / \
        f"{year}_estimation_GTI_{f}min_{method}_{pyrano_var_name}pyrano.csv"
    
    ###############################################################################
    
    # Load and filter data
    glob_estim_data = pd.read_csv(input_filename,parse_dates=['Timestamp'],index_col='Timestamp',
                             sep='\t',header=10)
    glob_estim_data = glob_estim_data.sort_index()
    glob_estim_data = glob_estim_data.loc[start_date:end_date]
    
    # Find the pyranometers removed for the beam/diffuse estimations
    pyrano_used = fct.find_content_between_braces(input_filename, line_number=5)
    pyrano_missing = fct.find_missing_elements(pyrano_used)
    
    ### Define the location
    # Read the first few lines of the file to find the location
    with open(input_filename, 'r') as file:
        header_lines = [next(file) for _ in range(6)]  # Read first 6 lines
    # Join the header lines and search for the location pattern
    header_text = ''.join(header_lines)
    match = re.search(r'Location:\s*(\d+\.\d+)\s*[NnSs]\s*(\d+\.\d+)\s*[EeWw]', header_text)
    lat_glob, lon_glob = float(match.groups()[0] ), float( match.groups()[1] )
    
    # Extract necessary columns
    beam = glob_estim_data['Beam'].values; diffuse = glob_estim_data['Diffuse'].values
    albedo = glob_estim_data['Albedo'].values
    beam_prime = glob_estim_data['Beam_prime']; diffuse_prime = glob_estim_data['Diffuse_prime']
    timestamps = glob_estim_data.index
    solar_angles = fct.calculate_solar_angles(timestamps, lat_glob, lon_glob)
    theta_z = solar_angles['zenith'].values; theta_z[theta_z > 88] = np.nan
    
    # Calculate irradiance for each combination of tilt and azimuth angles
    df_results = fct.calculate_GTI_for_orientations(
        solar_angles, tilt_angles_calc, azimuth_angles_calc, azimuth_angles,
        beam_prime, diffuse_prime, albedo, theta_z, lat_glob, lon_glob
    )
        
    # Set the index to the Timestamp column, filter out negative values and round the values
    df_results.set_index(timestamps, inplace=True)
    df_results = df_results.where(df_results >= 0, np.nan)
    df_results = df_results.round().astype('Int64')
    
    # Write the header and units to the file
    # Get the current date
    today = datetime.now().strftime("%Y-%m-%d")
    current_date = datetime.now().strftime("%Y-%m-%d")
    
    header = \
f"\
# Date of production: {current_date}\n\
# Produced by: Arthur Garreau\n\
# Estimation of solar irradiance for multiple orientations\
 (refer to Eq. 2.1 in Faiman et al. (1992)). Each column represents an\n\
# estimation for a specific tilt and azimuth, denoted as gti.azimuth_tilt.\n\
# Location: {lat_glob}N {lon_glob}E\n\
# The orientations missing are: {pyrano_missing}.\n\
# The orientations used are : {pyrano_used}\n\
#\n\
#\n\
# [UTC]\t[W m-2]\n"
    
    
    GTI_DATA_PATH.mkdir(parents=True, exist_ok=True)
    # Open the file and write the header and units
    with open(output_file, 'w', encoding='utf-8') as file:
        file.write(header)
    
    try:
        df_results.to_csv(output_file, mode= 'a', sep = "\t", na_rep ="NaN")
        print(f"Results saved to\n {output_file}")
    except Exception as e:
        print(f"Error saving results for\n {output_file}: {e}")


# %% CALC GTI with beam and diffuse from Ny-Ålesund
import pandas as pd
import numpy as np
from datetime import datetime
import glob_functions_calculation as fct
from config_path import NYA_DATA_PATH, GTI_DATA_PATH

# Parameters
year = 2025; f = 10 #min data frequency
############################## File Paths #####################################

bsrn_datafile = NYA_DATA_PATH / "NYA_radiation_2025-all.tab"

output_file = GTI_DATA_PATH / f"{year}_estimation_GTI_{f}min_NYA.csv"

###############################################################################

# Load the BSRN data
bsrn_data_full = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, 
                             parse_dates=['Date/Time'], index_col='Date/Time')
bsrn_data_full.index = bsrn_data_full.index.tz_localize('UTC')
bsrn_data_full = bsrn_data_full.rename_axis('Timestamp')
bsrn_data_full = bsrn_data_full.resample(f'{f}min').first()
if year == 2025: start_date, end_date = (f'{year}-03-16', f'{year}-07-31') 
bsrn_data_full = bsrn_data_full.loc[start_date:end_date]

lat_nya = 78.922700; lon_nya = 11.927300

# Extract necessary columns
beam = bsrn_data_full['DIR']; diffuse = bsrn_data_full['DIF']
albedo = bsrn_data_full['SWU']/bsrn_data_full['SWD']
filtered_albedo = albedo.where((albedo.index.hour >= 10) & (albedo.index.hour < 14))
daily_mean_albedo = filtered_albedo.resample('1D').mean()
new_albedo = albedo*np.nan
for day in daily_mean_albedo.index:
    daily_mean = daily_mean_albedo[str(day)]
    new_albedo[new_albedo.index.date == day.date()] = daily_mean
    
timestamps = bsrn_data_full.index
solar_angles = fct.calculate_solar_angles(timestamps, lat_nya, lon_nya)
theta_z = solar_angles['zenith']; theta_z[theta_z > 80] = np.nan
I_0 = pvlib.irradiance.get_extra_radiation(timestamps).values
diffuse_prime, beam_prime = fct.calculate_Dprime_and_Bprime(I_0, theta_z, diffuse, beam)

### Define the range of tilt and azimuth angles
tilt_angles_calc = np.arange(0, 181, 5)  
azimuth_angles_calc = np.arange(-180, 180, 15)
azimuth_angles_names = np.arange(0, 360, 15) 

# Calculate irradiance for each combination of tilt and azimuth angles
df_results = fct.calculate_GTI_for_orientations(
    solar_angles, tilt_angles_calc, azimuth_angles_calc, azimuth_angles,
    beam_prime, diffuse_prime, albedo, theta_z, lat_glob, lon_glob
)
# Set the index to the Timestamp column, filter out negative values and round the values
df_results.set_index(timestamps, inplace=True)
df_results = df_results.where(df_results >= 0, np.nan)
df_results = df_results.replace([np.inf, -np.inf], np.nan).round().astype('Int64')

today = datetime.now().strftime("%Y-%m-%d")
current_date = datetime.now().strftime("%Y-%m-%d")

header = \
f"\
# Date of production: {current_date}\n\
# Produced by: Arthur Garreau\n\
# Estimation of solar irradiance for multiple orientations \
(refer to Eq. 2.1 in Faiman et al. (1992)) with the BSRN beam and diffuse \
irradiance data of  Ny-Ålesund. Each column represents an \
estimation for a specific tilt and azimuth, denoted as gti.azimuth_tilt.\n\
# Location: {lat_nya}N {lon_nya}E\n\
#\n\
#\n\
#\n\
#\n\
#\n\
# [UTC]\t[W m-2]\n"

# Create the folder if it doesn't exist already
GTI_DATA_PATH.mkdir(parents=True, exist_ok=True)

# Open the file and write the header and units
with open(output_file, 'w', encoding='utf-8') as file:
    file.write(header)
try:
    df_results.to_csv(output_file, mode= 'a', sep = "\t", na_rep ="NaN")
    print(f"Results saved to\n {output_file}")
except Exception as e:
    print(f"Error saving results for\n {output_file}: {e}")
