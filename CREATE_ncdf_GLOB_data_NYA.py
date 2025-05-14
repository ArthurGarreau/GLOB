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

############################## File Paths #####################################

data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB")

###############################################################################

f = 10 #minute

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
        'declination': np.degrees(delta),
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


# ---- Data Preparation ---- #
file_2025 = data_path / "GLOB_data_30sec_2025_NYA.dat"

# Read and combine data
df_2025 = read_and_preprocess_data(file_2025)

combined_df = df_2025
combined_df.index = pd.to_datetime(combined_df.index)

# Resample to f-minute intervals and compute the mean
combined_df.index = combined_df.index + pd.Timedelta(minutes=f/2)
resampled_df = combined_df.resample(f'{f}min').mean()
combined_df.index = combined_df.index - pd.Timedelta(minutes=f/2)

# Convert DataFrame to Xarray Dataset
ds = xr.Dataset.from_dataframe(resampled_df)

timestamps_ds = pd.to_datetime(ds['Timestamp']).astype('datetime64[ns]')
timestamps_ds = timestamps_ds.tz_localize('UTC')

# Add metadata to the Timestamp variable
ds['Timestamp'] = timestamps_ds.values  # Important to add the .values to have the format YYYY-MM-DD hh:mm:ss
ds['Timestamp'].attrs.update({
    'calendar': 'gregorian',
    'long_name': 'UTC time',
    'standard_name': 'time'
})

ds = ds.drop_duplicates(dim='Timestamp')
timestamps_ds = timestamps_ds.drop_duplicates()

# Compute albedo and filter valid values (0 <= albedo <= 1)
albedo = ds['GHI_ground'] / ds['GHI']
albedo = albedo.where((albedo >= 0) & (albedo <= 1))

# Filter the albedo data between 10:00 and 12:00 for each day
filtered_albedo = albedo.sel(Timestamp=albedo.indexes['Timestamp'][albedo.indexes['Timestamp'].hour >= 10])
filtered_albedo = filtered_albedo.sel(Timestamp=filtered_albedo.indexes['Timestamp'][filtered_albedo.indexes['Timestamp'].hour < 12])

# Group by day and calculate the mean albedo for each day
daily_mean_albedo = filtered_albedo.resample(Timestamp='1D').mean(skipna=True)

# Create a new variable with the same time frequency as the original albedo data
new_albedo = xr.full_like(albedo, np.nan)  # Initialize with NaNs

# Set the daily mean albedo values for each day
for day in daily_mean_albedo.Timestamp.values:
    daily_mean = daily_mean_albedo.sel(Timestamp=day).item()
    new_albedo = new_albedo.where(new_albedo.Timestamp.dt.floor('D') != day, daily_mean)


# Add albedo to the original dataset
ds['albedo'] = (('Timestamp'), new_albedo.values)
ds['albedo'].attrs.update({
    'units': 'dimensionless',
    'long_name': 'Surface Albedo',
    'standard_name': 'albedo'
})

# ---- Solar Angles Calculation ---- #

# Define location coordinates
latitude = 78.92240
longitude = 11.92174

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

# Add latitude and longitude to the 5-minute dataset (ds)
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

output_file = data_path / f"GLOB_data_{f}min_2025.nc"

# Remove conflicting attributes from 'Timestamp' before saving
for attr in ['calendar', 'units']:
    if attr in ds['Timestamp'].attrs:
        del ds['Timestamp'].attrs[attr]

# Save datasets
ds.to_netcdf(output_file)

# Print variables in the dataset
for d in ds.variables:
    print(d)
print(f"{f}-minute NetCDF file created at: {output_file}")
