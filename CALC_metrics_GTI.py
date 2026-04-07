# -*- coding: utf-8 -*-
"""
Error Analysis of GTI (Global Tilted Irradiance) Estimations
===========================================================
This script evaluates the performance of GTI estimations using cross-validation with GLOB pyranometer data.
It calculates RMSE, nRMSE, MBE, and nMBE for each tilted direction and configuration,
and saves the results to Excel files for both GLOB and NY-Ålesund (NYA) estimations.

There is a second part in the script to calculate the metrics from the GTI estimations 
madw with the beam and diffuse measurement of the BSRN station in Ny-Ålesund.

Key Features:
-------------
- Loads GLOB pyranometer data and GTI estimation results.
- Calculates error metrics for each tilted direction and pyranometer configuration.
- Performs cross-validation by removing specific orientations and comparing estimations.
- Outputs results to Excel files for further analysis.

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: May 30, 2025
"""

import pandas as pd
import numpy as np
import xarray as xr
from config_path import DATA_PATH
import glob_functions_calculation as fct

# Parameters
f = 10 #min data frequency

glob_datafile = DATA_PATH / "GLOB_data" / "GLOB_data_10min_2025.nc"
# Load the data
glob_data = xr.open_dataset(glob_datafile, engine="h5netcdf").to_dataframe().tz_localize('UTC')
        

for method in ['linear', 'nonlinear']:
    ############################## File Paths #####################################
    output_file = DATA_PATH / "Estim_GTI" / f"2025_error_GTI_{method}.xlsx"
    ###############################################################################
  
    results_df = pd.DataFrame()
    for idx, pyr_nr in enumerate([3, 4, 5, 9, 13, 25]):
        pyrano_var = np.zeros(pyr_nr)
        GTI_estimation_datafile = DATA_PATH / "Estim_GTI" / \
            f"2025_estimation_GTI_{f}min_{method}_{pyr_nr}pyrano.csv"
        if pyr_nr == 3:
            if method == 'nonlinear':
                GTI_estimation_datafile = DATA_PATH / "Estim_GTI" / \
                    f"2025_estimation_GTI_{f}min_{method}_['GHI', 'N_90', 'N_45']pyrano.csv"
            if method == 'linear':
                GTI_estimation_datafile = DATA_PATH / "Estim_GTI" / \
                    f"2025_estimation_GTI_{f}min_{method}_['GHI', 'E_45', 'N_45']pyrano.csv"
        # Determine validation pyranometers
        removed_pyrano = fct.find_content_between_braces(GTI_estimation_datafile, line_number=4)
        validation_pyrano = [
            'NE_45', 'NE_90', 'NE_135', 'SE_45', 'SE_90', 'SE_135',
            'SW_45', 'SW_90', 'SW_135', 'NW_45', 'NW_90', 'NW_135'
        ]
    
        GTI_estimation_label = fct.create_gti_estimation_label(validation_pyrano)
        GTI_estimation_data = pd.read_csv(GTI_estimation_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)
    
        # Initialize a list to store error metrics
        error_metrics = []
    
        # Calculate metrics for each direction in validation_pyrano
        for direction in validation_pyrano:
            if direction in glob_data.columns and direction in GTI_estimation_label:
                GTI_glob = glob_data[direction]
                GTI_estimation_col = GTI_estimation_label[direction]
                GTI_estimation = GTI_estimation_data[GTI_estimation_col]
    
                GTI_glob, GTI_estimation = GTI_glob.align(GTI_estimation, join='inner', axis=0)
    
                nRMSE = int(fct.calculate_nrmse(GTI_glob, GTI_estimation))
                nMBE = int(fct.calculate_nmbe(GTI_glob, GTI_estimation))
                nMAE = int(fct.calculate_nmae(GTI_glob, GTI_estimation))  # Added nMAE calculation
                AVG = int(np.nanmean(GTI_glob[~np.isnan(GTI_estimation)]))

                if pyr_nr == 25:
                    error_metrics.append({
                        'Dir': direction,
                        'nRMSE'+"_GLOB" + str(pyr_nr): nRMSE,
                        'nMBE'+"_GLOB" + str(pyr_nr): nMBE,
                        'nMAE'+"_GLOB" + str(pyr_nr): nMAE,  # Added nMAE to the dictionary
                        'Avg': AVG
                    })
                else:
                    error_metrics.append({
                        'Dir': direction,
                        'nRMSE'+"_GLOB" + str(pyr_nr): nRMSE,
                        'nMBE'+"_GLOB" + str(pyr_nr): nMBE,
                        'nMAE'+"_GLOB" + str(pyr_nr): nMAE  # Added nMAE to the dictionary
                    })
    
        # Create a DataFrame to store the error metrics
        error_metrics_df = pd.DataFrame(error_metrics)
        error_metrics_df = error_metrics_df.set_index('Dir').transpose()
        
        results_df = pd.concat([results_df, error_metrics_df])
        
    
    # Save the error metrics to a CSV file
    with open(output_file, 'w') as file:
        file.write(f"Error GTI based on virtual sensing (cross-validation) removing the orientations : {removed_pyrano}.\n")
    
    results_df.to_excel(output_file, sheet_name='Sheet1', index=True)
    print(f"Error metrics saved to {output_file}")


# %% Metrics for NYÅ estimations of GTI

import pandas as pd
import numpy as np
from config_path import DATA_PATH
import glob_functions_calculation as fct

# Parameters
year = 2025; f = 10 #min data frequency

############################## File Paths #####################################
glob_datafile = DATA_PATH / "GLOB_data" / "GLOB_data_10min_2025.nc"
GTI_estimation_datafile = DATA_PATH / "Estim_GTI" / f"2025_estimation_GTI_{f}min_NYA_HayDavies.csv"

output_file = DATA_PATH / "Estim_GTI" / "2025_error_GTI_NYA_HayDavies.xlsx"
###############################################################################

results_df = pd.DataFrame()

# Determine validation pyranometers
removed_pyrano = fct.find_content_between_braces(GTI_estimation_datafile, line_number=4)
validation_pyrano = [
    'NE_45', 'NE_90', 'NE_135', 'SE_45', 'SE_90', 'SE_135',
    'SW_45', 'SW_90', 'SW_135', 'NW_45', 'NW_90', 'NW_135'
]

GTI_estimation_label = fct.create_gti_estimation_label(validation_pyrano)

# Load the data
glob_data =  fct.read_netcdf(glob_datafile).to_dataframe()
GTI_estimation_data = pd.read_csv(GTI_estimation_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)

# Initialize a list to store error metrics
error_metrics = []

# Calculate metrics for each direction in validation_pyrano
for direction in validation_pyrano:
    if direction in glob_data.columns and direction in GTI_estimation_label:
        GTI_glob = glob_data[direction].tz_localize('UTC')
        GTI_estimation_col = GTI_estimation_label[direction]
        GTI_estimation = GTI_estimation_data[GTI_estimation_col]

        GTI_glob, GTI_estimation = GTI_glob.align(GTI_estimation, join='inner', axis=0)

        RMSE = int(fct.calculate_rmse(GTI_glob, GTI_estimation))
        nRMSE = int(fct.calculate_nrmse(GTI_glob, GTI_estimation))
        MBE = int(fct.calculate_mbe(GTI_glob, GTI_estimation))
        nMBE = int(fct.calculate_nmbe(GTI_glob, GTI_estimation))
        AVG = int(np.nanmean(GTI_glob))

        error_metrics.append({
            'Dir': direction,
            'nRMSE': nRMSE,
            'nMBE': nMBE,
            'Avg': AVG
        })

# Create a DataFrame to store the error metrics
error_metrics_df = pd.DataFrame(error_metrics)
error_metrics_df = error_metrics_df.set_index('Dir').transpose()

results_df = pd.concat([results_df, error_metrics_df])

# Save the error metrics to a CSV file
with open(output_file, 'w') as file:
    file.write(f"Error GTI based on virtual sensing (cross-validation) removing the orientations : {removed_pyrano}.\n")

results_df.to_excel(output_file, sheet_name='Sheet1', index=True)
print(f"Error metrics saved to {output_file}")

