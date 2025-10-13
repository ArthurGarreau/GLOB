# -*- coding: utf-8 -*-
"""
Error Analysis of Beam and Diffuse Irradiance Estimations
=========================================================
This script evaluates the performance of different models for estimating beam and diffuse irradiance
using BSRN reference data from Ny-Ålesund, Svalbard. It calculates RMSE, nRMSE, MBE, and nMBE
for each model and configuration, and saves the results to an Excel file.

Key Features:
-------------
- Loads BSRN reference data and GLOB pyranometer estimations.
- Calculates beam and diffuse irradiance using Erbs, Perez, and Orgill-Holland models.
- Computes error metrics (RMSE, nRMSE, MBE, nMBE) for each model and GLOB configuration.
- Handles albedo data and calculates average albedo and MBE for nonlinear methods.
- Outputs results to an Excel file for further analysis.

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: May 30, 2025
"""

import pandas as pd
import numpy as np
from config_path import DATA_PATH
import pvlib

# Parameters
# method = 'linear'; year = 2025; f = 10 #min data frequency
three_pyranometers = 'no'
############################## File Paths #####################################
bsrn_datafile = DATA_PATH / "NYA_BSRN_data"  / "NYA_radiation_2025-all.tab"
    
if three_pyranometers == 'yes':
    # Case with 3 pyranometers
    output_file = DATA_PATH / "Estim_Beam_Diffuse" / f"{year}_error_beam_diffuse_{method}_3pyrano.xlsx"
    pyrano_vars = [
        ['GHI', 'N_90', 'N_45'],
        ['GHI', 'E_45', 'W_45'],
        ['GHI', 'S_45', 'W_45'],
        ['GHI', 'E_45', 'W_45']] 
else:
    pyrano_vars = [5, 9, 13, 25] # General case
    output_file = DATA_PATH / "Estim_Beam_Diffuse"  / f"{year}_error_beam_diffuse_{method}.xlsx"
###############################################################################

# Metric calculation functions
def calculate_rmse(y_true, y_pred):
    return np.sqrt(np.nanmean((y_true - y_pred) ** 2))

def calculate_nrmse(y_true, y_pred):
    rmse = calculate_rmse(y_true, y_pred)
    return rmse / (np.nanmean(y_true)) *100

def calculate_mbe(y_true, y_pred):
    return np.nanmean(y_pred - y_true)

def calculate_nmbe(y_true, y_pred):
    mbe = calculate_mbe(y_true, y_pred)
    return mbe / np.nanmean(y_true) *100

# Calculate metrics for each model    
def calculate_metrics(beam_est, diffuse_est, method_name):
    rmse_beam = int(calculate_rmse(merged_data['Beam_bsrn'], beam_est))
    nrmse_beam = int(calculate_nrmse(merged_data['Beam_bsrn'], beam_est))
    mbe_beam = int(calculate_mbe(merged_data['Beam_bsrn'], beam_est))
    nmbe_beam = int(calculate_nmbe(merged_data['Beam_bsrn'], beam_est))
    rmse_diffuse = int(calculate_rmse(merged_data['Diffuse_bsrn'], diffuse_est))
    nrmse_diffuse = int(calculate_nrmse(merged_data['Diffuse_bsrn'], diffuse_est))
    mbe_diffuse = int(calculate_mbe(merged_data['Diffuse_bsrn'], diffuse_est))
    nmbe_diffuse = int(calculate_nmbe(merged_data['Diffuse_bsrn'], diffuse_est))
    data = {
        'Method': method_name + str(pyr_nr),
        'Beam_RMSE': rmse_beam,
        'Beam_nRMSE': nrmse_beam,
        'Beam_MBE': mbe_beam,
        'Beam_nMBE': nmbe_beam,
        'Diffuse_RMSE': rmse_diffuse,
        'Diffuse_nRMSE': nrmse_diffuse,
        'Diffuse_MBE': mbe_diffuse,
        'Diffuse_nMBE': nmbe_diffuse,
        'Albedo_Avg': np.nan,
        'Albedo_MBE': np.nan
    }

    # Return a DataFrame with a single row
    return pd.DataFrame([data])

# Define the location
lat_glob, lon_glob = 78.9224, 11.92174


# Load and prepare data
bsrn_data = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col=['Date/Time'])
bsrn_data.index = bsrn_data.index.tz_localize('UTC')
albedo = bsrn_data['SWU'] / bsrn_data['SWD']
albedo = albedo.where((albedo.index.hour >= 10) & (albedo.index.hour < 14))
daily_mean_albedo = albedo.resample('1D').mean()
new_albedo = albedo*np.nan
for day in daily_mean_albedo.index:
    daily_mean = daily_mean_albedo[str(day)]
    new_albedo[new_albedo.index.date == day.date()] = daily_mean
    
timestamps = bsrn_data.index
bsrn_data['Albedo'] = new_albedo.loc[timestamps]
bsrn_data = bsrn_data.rename(columns={'DIR': 'Beam', 'DIF': 'Diffuse'})

# Calculate estimations with different models
ghi_bsrn = bsrn_data['SWD']
time = bsrn_data.index
solar_position = pvlib.solarposition.get_solarposition(time, lat_glob, lon_glob)
zenith = solar_position['zenith'].values

# Erbs model
erbs = pvlib.irradiance.erbs(ghi_bsrn, zenith, time)
beam_erbs, diffuse_erbs = erbs['dni'], erbs['dhi']

# Perez model
perez_dni = pvlib.irradiance.dirint(ghi_bsrn, zenith, time)
beam_perez, diffuse_perez = perez_dni, ghi_bsrn - np.cos(np.radians(zenith)) * perez_dni

# Orgill and Holland model
orgill_hollands = pvlib.irradiance.orgill_hollands(ghi_bsrn, zenith, time)
beam_oh, diffuse_oh = orgill_hollands['dni'], orgill_hollands['dhi']


results_df = pd.DataFrame()
for idx, pyr_nr in enumerate(pyrano_vars):

    ######################### Input File Paths ################################
    if three_pyranometers == 'yes':
        # Case with 3 pyranometers
        beam_diff_estim_datafile = DATA_PATH / "Estim_Beam_Diffuse" / \
            f"{year}_error_beam_diffuse_{method}_{pyr_nr}.xlsx"
    else: 
        beam_diff_estim_datafile = DATA_PATH / "Estim_Beam_Diffuse" / \
            f"{year}_estimation_beam_diffuse_{f}min_{method}_{str(pyr_nr)}pyrano.csv"
    ###########################################################################
  
    glob_estimation_data = pd.read_csv(beam_diff_estim_datafile, sep='\t', parse_dates=['Timestamp'], index_col=['Timestamp'], skiprows=10)

    # Align datasets
    merged_data = pd.merge(glob_estimation_data, bsrn_data, left_index=True, right_index=True, how='inner', suffixes=('_glob', '_bsrn'))
    columns_to_check = ['Beam_glob', 'Diffuse_glob', 'Beam_bsrn', 'Diffuse_bsrn']
    merged_data = merged_data.dropna(subset=columns_to_check)

  
    glob_df = calculate_metrics(merged_data['Beam_glob'], merged_data['Diffuse_glob'], 'GLOB')

    glob_df['Albedo_Avg'] = np.round(np.nanmean(merged_data['Albedo_glob']), 2)

    if method == "nonlinear":
        merge_albedo = merged_data.dropna(subset=['Albedo_glob', 'Albedo_bsrn'])
        glob_df['Albedo_Avg'] = np.round(np.nanmean(merge_albedo['Albedo_glob']), 2)
        glob_df['Albedo_MBE'] = np.round(calculate_mbe(merge_albedo['Albedo_bsrn'], merge_albedo['Albedo_glob']), 2)

    # Combine the DataFrames
    results_df = pd.concat([results_df, glob_df], ignore_index=True)


erbs_df = calculate_metrics(beam_erbs, diffuse_erbs, 'Erbs')
perez_df = calculate_metrics(beam_perez, diffuse_perez, 'Perez')
oh_df = calculate_metrics(beam_oh, diffuse_oh, 'Org_Hol')

# Combine the DataFrames
results_df = pd.concat([results_df, erbs_df, perez_df, oh_df], ignore_index=True)

# Write the results to a CSV file
with open(output_file, 'w') as file:
    file.write("Error of the beam and diffuse components (RMSE and MBE) based "
"on Ny-Alesund BSRN data.\n")
results_df.to_excel(output_file, sheet_name='Sheet1', index=False)
print(f"Results saved to {output_file}")
