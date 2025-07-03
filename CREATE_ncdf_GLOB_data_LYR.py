# -*- coding: utf-8 -*-
"""
GLOB Data Processing and NetCDF Conversion Script
=================================================

This script processes GLOB data from a CSV file, calculates solar angles, computes surface albedo,
and converts the data into a NetCDF file format.

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

# Add script path for importing the functions in glob_functions_Faiman.py
sys.path.append(str(SCRIPT_PATH))
import glob_functions_Faiman as fct

f = 10  # minute

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

def resample_data(df, frequency):
    """
    Resample data to a specified frequency and compute the mean.
    """
    df.index = df.index + pd.Timedelta(minutes=frequency / 2)
    resampled_df = df.resample(f'{frequency}min').mean()
    return resampled_df

def convert_to_xarray_with_metadata(resampled_df):
    """
    Convert a resampled DataFrame to an xarray Dataset and add metadata to the Timestamp variable.

    Parameters:
        - resampled_df: Pandas DataFrame that has been resampled.

    Returns:
        - xarray Dataset with added metadata for the Timestamp variable.
    """
    ds = xr.Dataset.from_dataframe(resampled_df)
    timestamps_ds = pd.to_datetime(ds['Timestamp']).astype('datetime64[ns]').tz_localize('UTC')

    ds['Timestamp'] = timestamps_ds.values
    ds['Timestamp'].attrs.update({
        'calendar': 'gregorian',
        'long_name': 'UTC time',
        'standard_name': 'time'
    })

    ds = ds.drop_duplicates(dim='Timestamp')
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
    
    timestamps_ds = pd.to_datetime(ds['Timestamp']).astype('datetime64[ns]').tz_localize('UTC')
    solar_angles = fct.calculate_solar_angles(timestamps_ds, latitude, longitude)

    for var, values in solar_angles.items():
        ds[var] = (('Timestamp'), values)
        ds[var].attrs.update({
            'units': 'degrees',
            'long_name': f'Solar {var.replace("_", " ")}',
            'standard_name': var
        })

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

def compute_and_filter_albedo(ds, ds_kZ):
    """
    Compute albedo from the irradiance of the K&Z in Adventdalen, filter valid values, and calculate daily mean albedo between 10:00 and 12:00.
    Data ds_KZ are from :
        Garreau, A., Shestov, A., Sikora, S., & Sjöblom, A. (2024). Aggregated SW and LW irradiance downwelling and upwelling in Adventdalen [NetCDF4-CF]. Arctic Data Centre. https://doi.org/10.21343/psy9-3e97

    Parameters:
        - ds: xarray Dataset containing the data.

    Returns:
        - xarray Dataset with added albedo data.
    """
    albedo = ds_kZ['SWup'] / ds_kZ['SWdown']
    albedo = albedo.drop_duplicates(dim='time')
    albedo = albedo.where((albedo >= 0) & (albedo <= 1))

    filtered_albedo = albedo.where((albedo.indexes['time'].hour >= 10) & (albedo.indexes['time'].hour < 12))
    daily_mean_albedo = filtered_albedo.resample(time='1D').mean(skipna=True)
    
    # The instrument was out of order in July 2024 and therefore the albedo is set to the value it had just before this event. 
    # The ground was free of snow during this period so the assumption is realistic.
    mask = (daily_mean_albedo['time'] >= pd.Timestamp('2024-07-08')) & (daily_mean_albedo['time'] <= pd.Timestamp('2024-08-19'))
    daily_mean_albedo = daily_mean_albedo.where(~mask, 0.115) # Assumed value for the missing data period.

    new_albedo = xr.full_like(ds['GHI'], np.nan)
    for day in daily_mean_albedo.time.values:
        daily_mean = daily_mean_albedo.sel(time=day).item()
        new_albedo = new_albedo.where(new_albedo.Timestamp.dt.floor('D') != day, daily_mean)

    aligned_albedo = new_albedo.interp(Timestamp=ds['Timestamp'].values, method='linear')
    ds['albedo'] = (('Timestamp'), aligned_albedo.values)
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




# %% ---- Process Data for each year ---- #
latitude = 78.200318; longitude = 15.840308 # Longyearbyen
file_KZ = DATA_PATH.parent / "Irradiance_ncdf" / "Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"

for year in [2023, 2024]:
    
    file_path = DATA_PATH / f"GLOB_data_30sec_{year}.dat"
    
    df = read_and_preprocess_data(file_path)
    
    df.index = pd.to_datetime(df.index)
    
    # Resample to f-minute intervals and compute the mean centered on the averaging period
    df.index = df.index + pd.Timedelta(minutes=f/2)
    resampled_df = df.resample(f'{f}min').mean()

    ds = convert_to_xarray_with_metadata(resampled_df)
    
    ds = add_solar_angles_and_coordinates(ds, fct, latitude, longitude)
    
    ds_kZ = xr.open_dataset(file_KZ)
    ds = compute_and_filter_albedo(ds, ds_kZ)
    
    ds = spectral_correction(ds)
    
    output_file = DATA_PATH / f"GLOB_data_{f}min_{year}.nc"
    save_to_netcdf(ds, output_file)
    ds.close(); ds_kZ.close()
