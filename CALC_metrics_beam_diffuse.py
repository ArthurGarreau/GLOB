# -*- coding: utf-8 -*-
"""
Created on Fri May 30 17:15:59 2025

@author: arthurg
"""
import pandas as pd
import numpy as np
from config_path import SCRIPT_PATH, DATA_PATH
import pvlib, re

# Parameters
# method = 'linear'; year = 2025; f = 10 #min data frequency

for i in [5, 9, 13, 25]:
    pyrano_var = np.zeros(i)
    ############################## File Paths #####################################
    bsrn_datafile = DATA_PATH.parent / "Irradiance" / "NYA" / "NYA_radiation_2025-all.tab"
    beam_diff_estim_datafile = SCRIPT_PATH / 'Data' / "Beam_Diffuse_Estimations" / \
        f"{year}_estimation_beam_diffuse_{f}min_{method}_{len(pyrano_var)}pyrano.csv"
    
    output_file = SCRIPT_PATH / 'Data' / "Beam_Diffuse_Estimations" / \
        f"{year}_error_beam_diffuse_{method}_{len(pyrano_var)}pyrano.csv"
    
    ###############################################################################
    
    # Define the location
    with open(beam_diff_estim_datafile, 'r') as file:
        header_lines = [next(file) for _ in range(10)]
    header_text = ''.join(header_lines)
    match = re.search(r'Location:\s*(\d+\.\d+)\s*[NnSs]\s*(\d+\.\d+)\s*[EeWw]', header_text)
    lat_glob, lon_glob = float(match.groups()[0]), float(match.groups()[1])
    
    # Metric calculation functions
    def calculate_rmse(y_true, y_pred):
        return np.sqrt(np.nanmean((y_true - y_pred) ** 2))
    
    def calculate_mbe(y_true, y_pred):
        return np.nanmean(y_true - y_pred)
    
    def calculate_mape(y_true, y_pred):
        non_zero_mask = (y_pred != 0)
        if non_zero_mask.any():
            return np.nanmean(np.abs((y_true[non_zero_mask] - y_pred[non_zero_mask]) / y_pred[non_zero_mask]) * 100)
        else:
            return np.nan
    
    # Load and prepare data
    bsrn_data = pd.read_csv(bsrn_datafile, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col=['Date/Time'])
    bsrn_data.index = bsrn_data.index.tz_localize('UTC')
    bsrn_data['Albedo'] = bsrn_data['SWU'] / bsrn_data['SWD']
    bsrn_data = bsrn_data.rename(columns={'DIR': 'Beam', 'DIF': 'Diffuse'})
    
    glob_estimation_data = pd.read_csv(beam_diff_estim_datafile, sep='\t', parse_dates=['Timestamp'], index_col=['Timestamp'], skiprows=10)
    
    # Align datasets
    merged_data = pd.merge(glob_estimation_data, bsrn_data, left_index=True, right_index=True, how='inner', suffixes=('_glob', '_bsrn'))
    columns_to_check = ['Beam_glob', 'Diffuse_glob', 'Beam_bsrn', 'Diffuse_bsrn']
    merged_data = merged_data.dropna(subset=columns_to_check)
    
    # Calculate estimations with different models
    ghi_bsrn = merged_data['SWD']
    time = merged_data.index
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
    
    # Calculate metrics for each model
    def calculate_metrics(beam_est, diffuse_est, method_name):
        avg_beam = int(np.nanmean(beam_est))
        avg_diffuse = int(np.nanmean(diffuse_est))
        rmse_beam = int(calculate_rmse(merged_data['Beam_bsrn'], beam_est))
        mbe_beam = int(calculate_mbe(merged_data['Beam_bsrn'], beam_est))
        mape_beam = int(calculate_mape(merged_data['Beam_bsrn'], beam_est))
        rmse_diffuse = int(calculate_rmse(merged_data['Diffuse_bsrn'], diffuse_est))
        mbe_diffuse = int(calculate_mbe(merged_data['Diffuse_bsrn'], diffuse_est))
        mape_diffuse = int(calculate_mape(merged_data['Diffuse_bsrn'], diffuse_est))
    
        data = {
            'Method': method_name,
            'Beam_Avg': avg_beam,
            'Beam_RMSE': rmse_beam,
            'Beam_MBE': mbe_beam,
            'Beam_MAPE': mape_beam,
            'Diffuse_Avg': avg_diffuse,
            'Diffuse_RMSE': rmse_diffuse,
            'Diffuse_MBE': mbe_diffuse,
            'Diffuse_MAPE': mape_diffuse,
            'Albedo_Avg':np.nan,
            'Albedo_MBE':np.nan
        }
        
        # Return a DataFrame with a single row
        return pd.DataFrame([data])
            
    glob_df = calculate_metrics(merged_data['Beam_glob'], merged_data['Diffuse_glob'], 'GLOB')
    
    glob_df['Albedo_Avg'] = np.round(np.nanmean(merged_data['Albedo_glob']),2)
    
    if method == "nonlinear":
        merge_albedo = merged_data.dropna(subset=['Albedo_glob', 'Albedo_bsrn'])
        glob_df['Albedo_Avg'] = np.round(np.nanmean(merge_albedo['Albedo_glob']),2),
        glob_df['Albedo_MBE'] = np.round(calculate_mbe(merge_albedo['Albedo_bsrn'], merge_albedo['Albedo_glob']),2),
    
    erbs_df = calculate_metrics(beam_erbs, diffuse_erbs, 'Erbs')
    perez_df = calculate_metrics(beam_perez, diffuse_perez, 'Perez')
    oh_df = calculate_metrics(beam_oh, diffuse_oh, 'Org_Hol')
    
    # Combine the DataFrames
    results_df = pd.concat([glob_df, perez_df, erbs_df, oh_df], ignore_index=True)
 
        
    # Write the results to a CSV file
    with open(output_file, 'w') as file:
        file.write("Error of the beam and diffuse components (RMSE and MBE) based on Ny-Alesund BSRN data.\n")
    results_df.to_csv(output_file, mode='a', sep='\t', index=False, header=True)
    print(f"Results saved to {output_file}")