# -*- coding: utf-8 -*-
"""
Parallel Beam and Diffuse Irradiance Calculation for GLOB
=========================================================

This script calculates beam and diffuse irradiance using the Faiman et al. (1992) method from GLOB measurements,
parallelizing the loop over pyrano_vars for improved performance.

Key Features:
-------------
- Loads GLOB data from a NetCDF file.
- Calculates solar angles and irradiance components using custom functions.
- Estimates beam and diffuse irradiance for each timestamp, in parallel.
- Saves results for each pyrano_var configuration.

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: February 18, 2026
"""

# %% Load Libraries
import pandas as pd
import xarray as xr
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from config_path import DATA_PATH
import glob_functions_calculation as fct


year = '2025'; f = 10  # min data frequency

############################## File Paths #####################################
glob_datafile = DATA_PATH / "GLOB_data" / f"GLOB_data_{f}min_{year}.nc"
output_file_path = DATA_PATH / "Estim_Beam_Diffuse"
###############################################################################

methods = ['linear', 'nonlinear']
pyrano_vars = [
    # ['GHI', 'N_45', 'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45', 'E_90', 'E_135',
    #  'SE_45', 'SE_90', 'SE_135', 'S_45', 'S_90', 'S_135', 'SW_45', 'SW_90', 'SW_135',
    #  'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135'],
    # ['GHI', 'N_45', 'N_90', 'N_135', 'E_45', 'E_90', 'E_135',
    #  'S_45', 'S_90', 'S_135', 'W_45', 'W_90', 'W_135'],
    # ['GHI', 'N_45', 'N_90', 'E_45', 'E_90', 'S_45', 'S_90', 'W_45', 'W_90'],
    # ['N_45', 'E_45', 'S_45', 'W_45'],
    ['GHI', 'E_45', 'S_45', 'W_45'],
    # ['GHI', 'N_90', 'N_45'],
    # ['GHI', 'E_45', 'N_45'],
    # ['GHI', 'S_45', 'W_45'],
]

def process_pyrano_var(pyrano_var, ds_glob, lat_glob, lon_glob, year, f, method, dates, output_file_path):
    """Process a single pyrano_var configuration and save results."""
    if len(pyrano_var) > 3:
        output_file = output_file_path / f"{year}_estimation_beam_diffuse_{f}min_{method}_{len(pyrano_var)}pyrano.csv"
    else:
        output_file = output_file_path / f"{year}_estimation_beam_diffuse_{f}min_{method}_{str(pyrano_var)}pyrano.csv"

    all_results = []
    for date in dates:
        date_str = date.strftime('%Y-%m-%d')
        df_glob_one_day = ds_glob.sel(Timestamp=date_str).to_pandas()
        df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')

        TIMESTAMPS = pd.date_range(
            start=f"{date_str} 00:00:00", end=f"{date_str} 23:50:00", freq='10min', tz='UTC')

        for timestamp in TIMESTAMPS:
            glob_value = df_glob_one_day.loc[timestamp]
            glob_value = glob_value.apply(pd.to_numeric, errors='coerce')
            solar_angles = fct.calculate_solar_angles(timestamp, lat_glob, lon_glob)
            new_result = fct.estimate_diffuse_beam(pyrano_var, glob_value, lat_glob, lon_glob, method=method)
            all_results.append([str(timestamp)] + new_result)
        print(date_str, pyrano_var)

    column_names = ['Timestamp', 'Diffuse', 'Beam', 'Albedo', 'Diffuse_prime', 'Beam_prime']
    df_all_results = pd.DataFrame(all_results, columns=column_names)
    df_all_results.set_index('Timestamp', inplace=True)

    current_date = datetime.now().strftime("%Y-%m-%d")
    header = (
        f"# Date of production: {current_date}\n"
        f"# Produced by: Arthur Garreau\n"
        f"# Estimation of beam and diffuse irradiance with GLOB using the least square error from the Faiman et al. (1992) method.\n"
        f"# In the linear method the albedo is measured, in the non linear the albedo is estimated.\n"
        f"# Location: {lat_glob}N {lon_glob}E\n"
        f"# The orientations used are : {pyrano_var}\n"
        f"#\n#\n#\n#[UTC]\t\t[W m-2]\t [W m-2]\t [/]\t [W m-2]\t [W m-2]\n"
    )

    with open(output_file, 'w', encoding='utf-8') as file:
        file.write(header)

    df_all_results.to_csv(output_file, index=True, mode='a', sep='\t', na_rep='NaN', encoding='utf-8')

def main():
    # Load GLOB data
    ds_glob = xr.open_dataset(glob_datafile, engine="h5netcdf")
    lat_glob = float(ds_glob.lat.values)
    lon_glob = float(ds_glob.lon.values)

    # Create a daily date range for the specified month and year
    if year == '2025':
        start_date, end_date = (f'{year}-03-16', f'{year}-09-05')
        dates = pd.date_range(start=start_date, end=end_date, freq='D')

    else:
        start_date_3, end_date_3 = ('2023-03-14', '2023-10-13')
        start_date_4, end_date_4 = ('2024-04-14', '2024-10-13')
        range_2023 = pd.date_range(start=start_date_3, end=end_date_3, freq='D')
        range_2024 = pd.date_range(start=start_date_4, end=end_date_4, freq='D')

        # Combine them into one index
        dates = range_2023.union(range_2024)

    # Parallel execution
    # --- Loop over Methods ---
    for method in methods:
        print(f"--- Starting processing for method: {method} ---")
        
        with ProcessPoolExecutor() as executor:
            executor.map(
                process_pyrano_var,
                pyrano_vars,
                [ds_glob] * len(pyrano_vars),
                [lat_glob] * len(pyrano_vars),
                [lon_glob] * len(pyrano_vars),
                [year] * len(pyrano_vars),
                [f] * len(pyrano_vars),
                [method] * len(pyrano_vars), # Passes current 'method' to all workers
                [dates] * len(pyrano_vars),
                [output_file_path] * len(pyrano_vars)
            )

    ds_glob.close()
if __name__ == "__main__":
    main()
