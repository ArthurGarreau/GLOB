# -*- coding: utf-8 -*-
"""
Error Analysis of Beam and Diffuse Irradiance Estimations
=========================================================
This script evaluates the performance of different models for estimating beam and diffuse irradiance
using BSRN reference data from Ny-Ålesund, Svalbard. It calculates nRMSE and nMBE
for each model and configuration, and saves the results to an Excel file.

Key Features:
-------------
- Loads BSRN reference data and GLOB pyranometer estimations.
- Calculates beam and diffuse irradiance using Erbs, Perez, Orgill-Holland, Muneer, Skartveit, and Spencer models.
- Computes error metrics (nRMSE, nMBE) for each model and GLOB configuration.
- Handles albedo data and calculates average albedo and MBE for nonlinear methods.
- Outputs results to an Excel file for further analysis.

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: May 30, 2025
"""

import pvlib
import pvmodels as pvm # Library for decomposition models
import pandas as pd
import numpy as np
from config_path import DATA_PATH
import glob_functions_calculation as fct

# Calculate metrics for each model
def calculate_metrics(beam_est, diffuse_est, method_name):
    nrmse_beam = int(fct.calculate_nrmse(merged_data['Beam_bsrn'], beam_est))
    nmbe_beam = int(fct.calculate_nmbe(merged_data['Beam_bsrn'], beam_est))
    nrmse_diffuse = int(fct.calculate_nrmse(merged_data['Diffuse_bsrn'], diffuse_est))
    nmbe_diffuse = int(fct.calculate_nmbe(merged_data['Diffuse_bsrn'], diffuse_est))

    # Calculate average of Beam and Diffuse from bsrn_data
    avg_beam_bsrn = int( np.nanmean(merged_data['Beam_bsrn']) )
    avg_diffuse_bsrn = int( np.nanmean(merged_data['Diffuse_bsrn']) )

    data = {
        'Method': method_name + str(pyr_nr),
        'Beam_nRMSE': nrmse_beam,
        'Beam_nMBE': nmbe_beam,
        'Diffuse_nRMSE': nrmse_diffuse,
        'Diffuse_nMBE': nmbe_diffuse,
        'Avg_Beam_bsrn': avg_beam_bsrn,
        'Avg_Diffuse_bsrn': avg_diffuse_bsrn,
        'Albedo_Avg': np.nan,
        'Albedo_MBE': np.nan
     }
                  
    # Return a DataFrame with a single row
    return pd.DataFrame([data])

# Parameters
f = 10 #min data frequency

for three_pyranometers in ['yes', 'no']:
    for method in ['linear', 'nonlinear']:   
############################## Output File Paths #####################################
        bsrn_datafile = DATA_PATH / "NYA_BSRN_data"  / "NYA_radiation_2025-all.tab"
        
        if three_pyranometers == 'yes':
            # Case with 3 pyranometers
            output_file = DATA_PATH / "Estim_Beam_Diffuse" / f"2025_error_beam_diffuse_{method}_3pyrano.xlsx"
            pyrano_vars = [
                ['GHI', 'N_90', 'N_45'],
                ['GHI', 'E_45', 'N_45'],
                ['GHI', 'S_45', 'W_45'],
                ['GHI', 'E_45', 'W_45']]
        else:
            pyrano_vars = [5, 9, 13, 25] # General case
            output_file = DATA_PATH / "Estim_Beam_Diffuse"  / f"2025_error_beam_diffuse_{method}.xlsx"
###############################################################################
        
        # Load and prepare data
        bsrn_data = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col=['Date/Time']).resample('10min').mean()
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
                
        results_df = pd.DataFrame()
        for idx, pyr_nr in enumerate(pyrano_vars):
        
######################### Input File Paths ####################################
            if three_pyranometers == 'yes':
                # Case with 3 pyranometers
                beam_diff_estim_datafile = DATA_PATH / "Estim_Beam_Diffuse" / \
                    f"2025_estimation_beam_diffuse_{f}min_{method}_{pyr_nr}pyrano.csv"
            else:
                beam_diff_estim_datafile = DATA_PATH / "Estim_Beam_Diffuse" / \
                    f"2025_estimation_beam_diffuse_{f}min_{method}_{str(pyr_nr)}pyrano.csv"
###############################################################################
        
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
                glob_df['Albedo_MBE'] = np.round(fct.calculate_mbe(merge_albedo['Albedo_bsrn'], merge_albedo['Albedo_glob']), 2)
        
            # Combine the DataFrames
            results_df = pd.concat([results_df, glob_df], ignore_index=True)
            

        
        # Define the location
        lat_glob, lon_glob = 78.9224, 11.92174
        
        
        
        # Calculate estimations with different models
        ghi_bsrn = bsrn_data['SWD']
        time = bsrn_data.index
        solar_position = pvlib.solarposition.get_solarposition(time, lat_glob, lon_glob)
        zenith = solar_position['zenith'].values; zenith[zenith > 85] = np.nan
        
        # Calculate Kt
        I_0 = pvlib.irradiance.get_extra_radiation(time) * np.cos(np.radians(zenith))
        kt = ghi_bsrn / I_0
        
        # Erbs model
        kd_erbs = pvm.erbs(kt)
        dhi_erbs = kd_erbs * ghi_bsrn
        dni_erbs = (ghi_bsrn - dhi_erbs) / np.cos(np.radians(zenith))
        
        # Perez model
        perez_dni = pvlib.irradiance.dirint(ghi_bsrn, zenith, time)
        dni_perez, dhi_perez = perez_dni, ghi_bsrn - np.cos(np.radians(zenith)) * perez_dni
        
        # Orgill and Holland model
        kd_orghol = pvm.orgill_holland(kt)
        dhi_orghol = kd_orghol * ghi_bsrn
        dni_orghol = (ghi_bsrn - dhi_orghol) / np.cos(np.radians(zenith))
        
        # Muneer model
        kd_muneer = pvm.muneer2(kt)
        dhi_muneer = kd_muneer * ghi_bsrn
        dni_muneer = (ghi_bsrn - dhi_muneer) / np.cos(np.radians(zenith))
        
        # Skartveit model
        kd_skartveit = pvm.skartveit1(kt, solar_position['elevation'].values)
        dhi_skartveit = kd_skartveit * ghi_bsrn
        dni_skartveit = (ghi_bsrn - dhi_skartveit) / np.cos(np.radians(zenith))
        
        # Spencer model
        kd_spencer = pvm.spencer(kt, lat_glob)
        dhi_spencer = kd_spencer * ghi_bsrn
        dni_spencer = (ghi_bsrn - dhi_spencer) / np.cos(np.radians(zenith))
               
        # Calculate metrics for each model
        perez_df = calculate_metrics(dni_perez, dhi_perez, 'Perez')
        erbs_df = calculate_metrics(dni_erbs, dhi_erbs, 'Erbs')
        oh_df = calculate_metrics(dni_orghol, dhi_orghol, 'Org_Hol')
        muneer_df = calculate_metrics(dni_muneer, dhi_muneer, 'Muneer2')
        mondol2_df = calculate_metrics(dni_skartveit, dhi_skartveit, 'Mondol2')
        spencer_df = calculate_metrics(dni_spencer, dhi_spencer, 'Spencer')
        
        # Combine the DataFrames
        results_df = pd.concat([results_df, spencer_df, muneer_df, mondol2_df, perez_df, oh_df, erbs_df], ignore_index=True)
        
        # Write the results to a CSV file
        with open(output_file, 'w') as file:
            file.write("Error of the beam and diffuse components (nRMSE and nMBE) based "
        "on Ny-Alesund BSRN data.\n")
        results_df.to_excel(output_file, sheet_name='Sheet1', index=False)
        print(f"Results saved to {output_file}")
    
################################ END ##########################################
