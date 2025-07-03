# -*- coding: utf-8 -*-
"""
GLOB Data Processing and NetCDF Conversion Script
=================================================

This script processes GLOB data from a CSV file, calculates solar angles, computes surface albedo,
resamples the data to 5-minute intervals, and converts the data into a NetCDF file format.

Key Features:
-------------
- Reads and preprocesses GLOB data from a CSV file.
- Calculates solar angles (zenith, elevation, azimuth, hour angle, declination) using pvlib.
- Computes surface albedo and filters valid values.
- Resamples the data to 5-minute intervals.
- Adds latitude and longitude information to the NetCDF file.
- Saves the processed data as a NetCDF file with appropriate metadata.

Dependencies:
-------------
- pandas
- xarray
- numpy
- pathlib

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 1, 2024
"""

import pandas as pd
import xarray as xr
import numpy as np
import sys
from config_path import SCRIPT_PATH, DATA_PATH
# Add script path for the importing the functions in glob_functions_Faiman.py
sys.path.append(str(SCRIPT_PATH))
import glob_functions_calculation as fct

f = int(input("Please enter the data frequency required in minutes: ")) #minute

# %% ---- Function Definitions ---- #

def read_and_preprocess_data(file_path):
    """
    Read and preprocess GLOB data from CSV.
    """
    df = pd.read_csv(file_path, skiprows=3)
    df.replace(-9999, np.nan, inplace=True)

    # Convert 'Timestamp' column to datetime and set it as index
    df['Timestamp'] = pd.to_datetime(df['Timestamp'].values, unit='ns')
    df.set_index('Timestamp', inplace=True)

    return df

def convert_df_to_xarray_with_metadata(resampled_df):
    """
    Convert a resampled DataFrame to an xarray Dataset and add metadata to the Timestamp variable.

    Parameters:
        - resampled_df: Pandas DataFrame that has been resampled.

    Returns:
        - xarray Dataset with added metadata for the Timestamp variable.
    """
    # Convert DataFrame to Xarray Dataset
    ds = xr.Dataset.from_dataframe(resampled_df)

    # Process timestamps
    timestamps_ds = pd.to_datetime(ds['Timestamp']).astype('datetime64[ns]').tz_localize('UTC')

    # Add metadata to the Timestamp variable
    ds['Timestamp'] = timestamps_ds.values
    ds['Timestamp'].attrs.update({
        'calendar': 'gregorian',
        'long_name': 'UTC time',
        'standard_name': 'time'
    })

    # Drop duplicates
    ds = ds.drop_duplicates(dim='Timestamp')
    timestamps_ds = timestamps_ds.drop_duplicates()

    return ds

def add_solar_angles_and_coordinates(ds, fct, latitude, longitude):
    """
    Add latitude, longitude, and solar angles to an xarray Dataset.

    Parameters:
        - ds: xarray Dataset containing the data.
        - fct: Module containing the function to calculate solar angles.
        - latitude: Latitude of the location (default is 78.92240).
        - longitude: Longitude of the location (default is 11.92174).

    Returns:
        - xarray Dataset with added latitude, longitude, and solar angles.
    """
    # Calculate solar angles
    timestamps_ds = pd.to_datetime(ds['Timestamp']).astype('datetime64[ns]').tz_localize('UTC')
    solar_angles = fct.calculate_solar_angles(timestamps_ds, latitude, longitude)

    # Add solar angles to the dataset
    for var, values in solar_angles.items():
        ds[var] = (('Timestamp'), values)
        ds[var].attrs.update({
            'units': 'degrees',
            'long_name': f'Solar {var.replace("_", " ")}',
            'standard_name': var
        })

    # Add latitude and longitude to the dataset
    ds['latitude'] = ((), latitude)
    ds['longitude'] = ((), longitude)

    ds['latitude'].attrs.update({
        'units': 'degrees',
        'long_name': 'Latitude of the location in Ny-Alesund',
        'standard_name': 'latitude'
    })

    ds['longitude'].attrs.update({
        'units': 'degrees',
        'long_name': 'Longitude of the location in Ny-Alesund',
        'standard_name': 'longitude'
    })

    return ds

def compute_and_filter_albedo(ds):
    """
    Compute albedo from GHI_ground and GHI, filter valid values, and calculate daily mean albedo between 10:00 and 12:00.

    Parameters:
        - ds: xarray Dataset containing the data.

    Returns:
        - xarray Dataset with added albedo data.
    """
    # Compute albedo and filter valid values (0 <= albedo <= 1)
    albedo = ds['GHI_ground'] / ds['GHI']
    albedo = albedo.where((albedo >= 0) & (albedo <= 1))

    # Filter the albedo data between 10:00 and 12:00 for each day
    filtered_albedo = albedo.where((albedo.indexes['Timestamp'].hour >= 10) & (albedo.indexes['Timestamp'].hour < 12))

    # Group by day and calculate the mean albedo for each day
    daily_mean_albedo = filtered_albedo.resample(Timestamp='1D').mean(skipna=True)

    # Create a new variable with the same time frequency as the original albedo data
    new_albedo = xr.full_like(albedo, np.nan)  # Initialize with NaNs

    # Set the daily mean albedo values for each day
    for day in daily_mean_albedo.Timestamp:
        daily_mean = daily_mean_albedo.sel(Timestamp=day).item()
        new_albedo.loc[new_albedo.Timestamp.dt.floor('D') == day] = daily_mean

    # Add albedo to the original dataset
    ds['albedo'] = new_albedo
    ds['albedo'].attrs.update({
        'units': 'dimensionless',
        'long_name': 'Surface Albedo',
        'standard_name': 'albedo'
    })

    return ds

def spectral_correction(ds, h=0.01):
    """
    Apply spectral correction to the silicon cell pyranometer data as explained in Balthazar et al. (2015).

    Parameters:
        - ds: xarray Dataset containing the data.
        - h: Height above sea level in km (default is 0.01).

    Returns:
        - xarray Dataset with corrected components.
    """
    # Extract zenith angle from the dataset
    theta_z = ds['zenith']

    # Calculate absolute air mass (AM_a)
    AM_a = np.exp(-0.0001184 * h) * (np.cos(np.radians(theta_z)) + 0.5057 * (96.080 - theta_z)**(-1.634))**(-1)

    # Calculate f1
    f1 = 0.000263 * (AM_a)**3 - 0.00632 * (AM_a)**2 + 0.054 * (AM_a) + 0.932

    # Calculate f2
    f2 = -0.00000004504 * theta_z**3 - 0.00001357 * theta_z**2 + 0.0006074 * theta_z + 1

    # List of components to correct
    components = ['GHI', 'S_45', 'S_90', 'S_135', 'SW_45', 'SW_90', 'SW_135',
                  'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135', 'N_45',
                  'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45', 'E_90',
                  'E_135', 'SE_45', 'SE_90', 'SE_135']

    # Apply the correction to each component
    for component in components:
        if component in ds:
            ds[f'{component}'] = ds[component] / (f1 * f2)

    return ds

def save_to_netcdf(ds, output_file):
    """
    Save the xarray Dataset to a NetCDF file.
    """
    for attr in ['calendar', 'units']:
        if attr in ds['Timestamp'].attrs:
            del ds['Timestamp'].attrs[attr]

    ds.to_netcdf(output_file)
    print(f"{f}-minute NetCDF file created at: {output_file}")

# %% ---- Process Data ---- #
file_2025 = DATA_PATH / "GLOB_data_30sec_2025_NYA.dat"

# Read and combine data
df = read_and_preprocess_data(file_2025)

df.index = pd.to_datetime(df.index)

# Resample to f-minute intervals and compute the mean centered on the averaging period
df.index = df.index + pd.Timedelta(minutes=f/2)
resampled_df = df.resample(f'{f}min').mean()
# resampled_df.index = resampled_df.index - pd.Timedelta(minutes=f/2)

ds = convert_df_to_xarray_with_metadata(resampled_df)

latitude=78.92240; longitude=11.92174 # Ny-Ålesund
ds = add_solar_angles_and_coordinates(ds, fct, latitude, longitude)

ds = compute_and_filter_albedo(ds)

ds = spectral_correction(ds)

output_file = DATA_PATH / f"GLOB_data_{f}min_2025.nc"
save_to_netcdf(ds, output_file)

# Print variables in the dataset
print('The ncdf dataset is from ', str(ds['Timestamp'].values[0])[0:10], 'to'
      , str(ds['Timestamp'].values[-1])[0:10])

