# -*- coding: utf-8 -*-
"""
Created on Fri May 30 17:15:59 2025

@author: arthurg
"""
import pandas as pd
import xarray as xr
import numpy as np
from config_path import SCRIPT_PATH, DATA_PATH
import glob_functions_calculation as fct

# method = 'nonlinear'; year = 2025; f = 10 #min data frequency

for i in [5, 9, 13, 25]:
    pyrano_var = np.zeros(i)
    ############################## File Paths #####################################
    glob_datafile = DATA_PATH / f"GLOB_data_10min_{year}.nc"
    gti_estimation_datafile = SCRIPT_PATH / 'Data' / "GTI_Estimations" / \
        f"{year}_estimation_GTI_{f}min_{method}_{len(pyrano_var)}pyrano.csv"
    
    output_file = SCRIPT_PATH / 'Data' / "GTI_Estimations" / \
        f"{year}_error_GTI_{method}_{len(pyrano_var)}pyrano.csv"
    
    ###############################################################################
    
    # Determine validation pyranometers
    removed_pyrano = fct.find_content_between_braces(gti_estimation_datafile, line_number=4)
    if removed_pyrano == ['']:
        validation_pyrano = ['S_90', 'N_90', 'E_90', 'W_90']
    else:
        validation_pyrano = [
            'NE_45', 'NE_90', 'NE_135', 'SE_45', 'SE_90', 'SE_135',
            'SW_45', 'SW_90', 'SW_135', 'NW_45', 'NW_90', 'NW_135'
        ]
    
    gti_estimation_label = fct.create_gti_estimation_label(validation_pyrano)
    
    # Load the data
    glob_data = xr.open_dataset(glob_datafile).to_dataframe()
    gti_estimation_data = pd.read_csv(gti_estimation_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)
    
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
    
    # Initialize a list to store error metrics
    error_metrics = []
    
    # Calculate metrics for each direction in validation_pyrano
    for direction in validation_pyrano:
        if direction in glob_data.columns and direction in gti_estimation_label:
            GTI_glob = glob_data[direction].tz_localize('UTC')
            gti_estimation_col = gti_estimation_label[direction]
            gti_estimation = gti_estimation_data[gti_estimation_col]
    
            GTI_glob, gti_estimation = GTI_glob.align(gti_estimation, join='inner', axis=0)
    
            RMSE = np.round(calculate_rmse(GTI_glob, gti_estimation))
            MBE = np.round(calculate_mbe(GTI_glob, gti_estimation))
            MAPE = np.round(calculate_mape(GTI_glob, gti_estimation))
            AVG = np.round(np.nanmean(GTI_glob))
    
            error_metrics.append({
                'Dir': direction,
                'Avg': AVG,
                'RMSE': RMSE,
                'MBE': MBE,
                'MAPE': MAPE,            
                'Units': 'W/m2'
            })
    
    # Create a DataFrame to store the error metrics
    error_metrics_df = pd.DataFrame(error_metrics)
    error_metrics_df = error_metrics_df.set_index('Dir').transpose()
    
    # Save the error metrics to a CSV file
    with open(output_file, 'w') as file:
        file.write(f"Error GTI based on virtual sensing (cross-validation) removing the orientations : {removed_pyrano}.\n")
    
    error_metrics_df.to_csv(output_file, mode='a', index=True, sep='\t')
    print(f"Error metrics saved to {output_file}")
