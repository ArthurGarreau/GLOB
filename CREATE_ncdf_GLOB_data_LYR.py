# -*- coding: utf-8 -*-
"""
GLOB Data Processing and NetCDF Conversion Script
=================================================

This script processes GLOB data from a CSV file, calculates solar angles, computes surface albedo,
and converts the data into a NetCDF file format. It includes functionalities for reading and preprocessing
data, calculating solar angles using pvlib, and adding metadata to the NetCDF file.

Key Features:
-------------
- Reads and preprocesses GLOB data from a CSV file.
- Calculates solar angles (zenith, elevation, azimuth, hour angle, declination) using pvlib.
- Computes surface albedo and filters valid values.
- Adds latitude and longitude information to the NetCDF file.
- Saves the processed data as a NetCDF file with appropriate metadata.

Dependencies:
-------------
- pandas
- xarray
- numpy
- pvlib
- pathlib

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 1, 2024
"""

import pandas as pd
import xarray as xr
import numpy as np
import pvlib
from pathlib import Path


# ---- Function Definitions ---- #

def calculate_solar_angles(timestamps, latitude, longitude, altitude=6, temperature=-6):
    """
    Recalculate solar angles for given timestamps.
    """
    solar_position = pvlib.solarposition.get_solarposition(
        time=timestamps, latitude=latitude, longitude=longitude, altitude=altitude, temperature=temperature
    )
    omega = np.radians(pvlib.solarposition.hour_angle(timestamps, longitude, solar_position.equation_of_time))
    delta = pvlib.solarposition.declination_spencer71(timestamps.day_of_year)

    angles = {
        'zenith_angle': solar_position['apparent_zenith'],
        'elevation_angle': solar_position['apparent_elevation'],
        'azimuth_angle': solar_position['azimuth'],
        'hour_angle': omega,
        'declination': delta,
    }
    return angles


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


# ---- File Paths ---- #

data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB")
file_2023 = data_path / "GLOB_data_30sec_2023.dat"
file_2024 = data_path / "GLOB_data_30sec_2024.dat"

file_KZ = data_path.parent / "Irradiance_ncdf" / "Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"

# %% ---- Create 30-sec dataset ---- #

# Read and combine data
df_2023 = read_and_preprocess_data(file_2023)
df_2024 = read_and_preprocess_data(file_2024)

combined_df = pd.concat([df_2023, df_2024])

# Convert DataFrames to Xarray Datasets
ds = xr.Dataset.from_dataframe(combined_df)

timestamps_ds = pd.to_datetime(ds['Timestamp']).astype('datetime64[ns]')
timestamps_ds = timestamps_ds.tz_localize('UTC')

# Add metadata to the Timestamp variable for both datasets

ds['Timestamp'] = timestamps_ds.values # Important to add the .values to have the format YYYY-MM-DD hh:mm:ss
ds['Timestamp'].attrs.update({
    'calendar': 'gregorian',
    'long_name': 'UTC time',
    'standard_name': 'time'
})

ds = ds.drop_duplicates(dim='Timestamp')
timestamps_ds = timestamps_ds.drop_duplicates()


# ---- Solar Angles Calculation ---- #

# Define location coordinates
latitude = 78.200318
longitude = 15.840308

solar_angles = calculate_solar_angles(timestamps_ds, latitude, longitude)

# Add solar angles to the dataset
for var, values in solar_angles.items():
    ds[var] = (('Timestamp'), values)
    ds[var].attrs.update({
        'units': 'degrees',
        'long_name': f'Solar {var.replace("_", " ")}',
        'standard_name': var
    })

# ---- Add Latitude and Longitude to NetCDF ---- #

# Add latitude and longitude to the 30-second dataset (ds)
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

# ---- Save to NetCDF ---- #

output_file_30sec = data_path / "GLOB_data_30sec_2023-24.nc"

# Remove conflicting attributes from 'Timestamp' before saving
for attr in ['calendar', 'units']:
    if attr in ds['Timestamp'].attrs:
        del ds['Timestamp'].attrs[attr]


# Save datasets
ds.to_netcdf(output_file_30sec)

print(f"30-second NetCDF file created at: {output_file_30sec}")


for d in ds.variables:
    print(d)

# %% ---- Create 5-min dataset ---- #

# Read and combine data
df_2023 = read_and_preprocess_data(file_2023)
df_2024 = read_and_preprocess_data(file_2024)

combined_df = pd.concat([df_2023, df_2024])
combined_df_5min = combined_df.resample('5min').mean()

ds_5min = xr.Dataset.from_dataframe(combined_df_5min)

timestamps_ds = pd.to_datetime(ds_5min['Timestamp']).astype('datetime64[ns]')
timestamps_ds = timestamps_ds.tz_localize('UTC')

# Add metadata to the Timestamp variable for both datasets

ds_5min['Timestamp'] = timestamps_ds.values # Important to add the .values to have the format YYYY-MM-DD hh:mm:ss
ds_5min['Timestamp'].attrs.update({
    'calendar': 'gregorian',
    'long_name': 'UTC time',
    'standard_name': 'time'
})

ds_5min = ds_5min.drop_duplicates(dim='Timestamp')
timestamps_ds = timestamps_ds.drop_duplicates()

# ---- Albedo Calculation ---- #

# Load KZ dataset
ds_kZ = xr.open_dataset(file_KZ)
# Compute albedo and filter valid values (0 <= albedo <= 1)
albedo = ds_kZ['SWup'] / ds_kZ['SWdown']
albedo = albedo.where(ds_kZ['SWup_quality'] == 'ok').drop_duplicates(dim='time')
albedo = albedo.where((albedo >= 0) & (albedo <= 1))


# Align albedo to match the 5-minute dataset's timestamps
aligned_albedo = albedo.interp(time=timestamps_ds.values, method='linear')
ds_5min['albedo'] = (('Timestamp'), aligned_albedo.values)
ds_5min['albedo'].attrs.update({
    'units': 'dimensionless',
    'long_name': 'Surface Albedo',
    'standard_name': 'albedo'
})


# ---- Solar Angles Calculation ---- #

# Define location coordinates
latitude = 78.200318
longitude = 15.840308

solar_angles = calculate_solar_angles(timestamps_ds, latitude, longitude)

# Add solar angles to the dataset
for var, values in solar_angles.items():
    ds_5min[var] = (('Timestamp'), values)
    ds_5min[var].attrs.update({
        'units': 'degrees',
        'long_name': f'Solar {var.replace("_", " ")}',
        'standard_name': var
    })

# ---- Add Latitude and Longitude to NetCDF ---- #

# Add latitude and longitude to the 30-second dataset (ds)
ds_5min['latitude'] = ((), latitude)
ds_5min['longitude'] = ((), longitude)

ds_5min['latitude'].attrs.update({
    'units': 'degrees',
    'long_name': 'Latitude of the location in Ny-Alesund',
    'standard_name': 'latitude'
})
ds_5min['longitude'].attrs.update({
    'units': 'degrees',
    'long_name': 'Longitude of the location in Ny-Alesund',
    'standard_name': 'longitude'
})

# ---- Save to NetCDF ---- #

output_file_5min = data_path / "GLOB_data_5min_2023-24.nc"

# Remove conflicting attributes from 'Timestamp' before saving
for attr in ['calendar', 'units']:
    if attr in ds_5min['Timestamp'].attrs:
        del ds_5min['Timestamp'].attrs[attr]

# Save datasets
ds_5min.to_netcdf(output_file_5min)

print(f"5-minute NetCDF file created at: {output_file_5min}")

for d in ds_5min.variables:
    print(d)

