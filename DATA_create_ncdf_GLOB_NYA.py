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
- Adds latitude and longitude information to the NetCDF file.
- Saves the processed data as a NetCDF file with appropriate metadata.

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 1, 2024
"""

import pandas as pd
import xarray as xr
import numpy as np
from config_path import DATA_PATH
import glob_functions_calculation as fct


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

def convert_df_to_xarray_with_metadata(resampled_df, latitude, longitude):
    """
    Convert a resampled DataFrame to an xarray Dataset and add metadata to the Timestamp variable.

    Parameters:
        - resampled_df: Pandas DataFrame that has been resampled.

    Returns:
        - xarray Dataset with added metadata for the Timestamp variable.
    """
    # Convert DataFrame to Xarray Dataset
    ds = xr.Dataset.from_dataframe(resampled_df)
    timestamps_ds = pd.to_datetime(ds['Timestamp']).astype('datetime64[ns]').tz_localize('UTC')
    ds['Timestamp'] = timestamps_ds.values
    ds['Timestamp'].attrs.update({
        'calendar': 'gregorian',
        'long_name': 'UTC time',
        'standard_name': 'time',
    })
    ds = ds.drop_duplicates(dim='Timestamp')
    ds['RECORD'] = ds['RECORD'].astype('int32')
    
    # Define metadata for each variable
    variable_metadata = {
        'RECORD': {
            'long_name': 'Record number',
            'units': '1',
        },
        'Temp': {
            'long_name': 'Temperature',
            'standard_name': 'air_temperature',
            'units': 'Celsius',
        },
        'GHI': {
            'long_name': 'Global horizontal irradiance',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'S_45': {
            'long_name': 'South-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'S_90': {
            'long_name': 'South-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'S_135': {
            'long_name': 'South-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'SW_45': {
            'long_name': 'Southwest-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'SW_90': {
            'long_name': 'Southwest-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'SW_135': {
            'long_name': 'Southwest-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'W_45': {
            'long_name': 'West-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'W_90': {
            'long_name': 'West-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'W_135': {
            'long_name': 'West-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'NW_45': {
            'long_name': 'Northwest-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'NW_90': {
            'long_name': 'Northwest-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'NW_135': {
            'long_name': 'Northwest-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'N_45': {
            'long_name': 'North-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'N_90': {
            'long_name': 'North-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'N_135': {
            'long_name': 'North-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'NE_45': {
            'long_name': 'Northeast-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'NE_90': {
            'long_name': 'Northeast-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'NE_135': {
            'long_name': 'Northeast-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'E_45': {
            'long_name': 'East-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'E_90': {
            'long_name': 'East-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'E_135': {
            'long_name': 'East-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'SE_45': {
            'long_name': 'Southeast-facing irradiance at 45° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'SE_90': {
            'long_name': 'Southeast-facing irradiance at 90° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'SE_135': {
            'long_name': 'Southeast-facing irradiance at 135° tilt',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'GHI_ground': {
            'long_name': 'Ground reflected horizontal irradiance',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
        },
        'albedo': {
            'long_name': 'Surface albedo',
            'standard_name': 'surface_albedo',
            'units': '1',
        },
    }
    
    # Add latitude and longitude as scalar coordinates
    ds.coords['lat'] = ((), latitude, {'units': 'degrees_north', 'long_name': 'Latitude of the location', 'standard_name': 'latitude'})
    ds.coords['lon'] = ((), longitude, {'units': 'degrees_east', 'long_name': 'Longitude of the location', 'standard_name': 'longitude'})

    # Assign metadata to each variable
    for var_name, attrs in variable_metadata.items():
        if var_name in ds:
            ds[var_name].attrs.update(attrs)

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
        ds[var] = (('Timestamp',), values.values)
        ds[var].attrs.update({
            'units': 'degrees',
            'long_name': f'Solar {var.replace("_", " ")}',
            'standard_name': 'solar' + var
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
    filtered_albedo = filtered_albedo.sortby('Timestamp')
    # Group by day and calculate the mean albedo for each day
    daily_mean_albedo = filtered_albedo.resample(Timestamp='1D').mean(skipna=True)

    # Create a new variable with the same time frequency as the original albedo data
    new_albedo = xr.full_like(albedo, np.nan)  # Initialize with NaNs

    # Set the daily mean albedo values for each day
    for day in daily_mean_albedo.Timestamp:
        daily_mean = daily_mean_albedo.sel(Timestamp=day).item()
        new_albedo.loc[dict(Timestamp=new_albedo.Timestamp.dt.floor('D') == day)] = daily_mean

    # Add albedo to the original dataset
    ds['albedo'] = (('Timestamp',), new_albedo.values)
    ds['albedo'].attrs.update({
        'units': '1',
        'long_name': 'Surface Albedo',
        'standard_name': 'surface_albedo'
    })

    return ds

   
def substitute_ghi_ground(df, df_bsrn):
    """
    Substitute GHI_ground values in df with SWU values from df_bsrn for the period 2025-05-09 to 2025-08-21.
    This is necessary because the GHI_ground sensor of GLOB fell off and therefore we have to replace those data
    by the BSRN data that presents ground reflected accurately measured next to GLOB.
    """
    # Ensure both dataframes have datetime indices
    df.index = pd.to_datetime(df.index)
    df_bsrn.index = pd.to_datetime(df_bsrn.index)

    # Extract the relevant period and column from df_bsrn
    # Assuming 'SWU' is the correct column for ground reflected data; replace if needed
    swu_values = df_bsrn.loc['2025-05-09':'2025-07-31', 'SWU']


    # Update df['GHI_ground'] with the extracted values
    df.loc[swu_values.index, 'GHI_ground'] = swu_values.values

    # Print some results for verification
    return df

def add_global_attributes(ds, **kwargs):
    """
    Add comprehensive global attributes to an xarray Dataset, following ACDD 1.3 conventions.

    Parameters:
        - ds: xarray Dataset to which attributes will be added.
        - **kwargs: Keyword arguments for all global attributes (e.g., title, summary, keywords, etc.).

    Returns:
        - xarray Dataset with added global attributes.
    """
    # Required ACDD attributes
    ds.attrs.update({
        'title': kwargs.get('title', ''),
        'summary': kwargs.get('summary', ''),
        'keywords': kwargs.get('keywords', ''),
        'Conventions': kwargs.get('Conventions', 'CF-1.8, ACDD-1.3'),
        'data_type': kwargs.get('data_type', 'netCDF-4'),
        'geospatial_lat_min': kwargs.get('geospatial_lat_min', ''),
        'geospatial_lat_max': kwargs.get('geospatial_lat_max', ''),
        'geospatial_lon_min': kwargs.get('geospatial_lon_min', ''),
        'geospatial_lon_max': kwargs.get('geospatial_lon_max', ''),
        'time_coverage_start': kwargs.get('time_coverage_start', ''),
        'time_coverage_end': kwargs.get('time_coverage_end', ''),
        'standard_name_vocabulary': kwargs.get('standard_name_vocabulary', 'CF Standard Name Table v82'),
        'featureType': kwargs.get('featureType', 'timeSeries'),
        'date_created': kwargs.get('date_created', ''),
        'history': kwargs.get('history', ''),
        'comments': kwargs.get('comments', ''),
        'creator_type': kwargs.get('creator_type', ''),
        'creator_institution': kwargs.get('creator_institution', ''),
        'creator_name': kwargs.get('creator_name', ''),
        'creator_url': kwargs.get('creator_url', ''),
        'creator_email': kwargs.get('creator_email', ''),
        'institution': kwargs.get('institution', ''),
        'license': kwargs.get('license', ''),
        'iso_topic_category': kwargs.get('iso_topic_category', ''),
        'publisher_name': kwargs.get('publisher_name', ''),
        'publisher_institution': kwargs.get('publisher_institution', ''),
        'publisher_url': kwargs.get('publisher_url', ''),
        'publisher_email': kwargs.get('publisher_email', ''),
        'publisher_type': kwargs.get('publisher_type', 'institution'),
        'station_name': kwargs.get('station_name', ''),
        'instrument_type': kwargs.get('instrument_type', ''),
    })

    return ds

def convert_int64_to_int32(ds):
    """
    Convert all int64 variables in the dataset to int32 for THREDDS compatibility.
    """
    for var in ds.data_vars:
        if ds[var].dtype == 'int64':
            ds[var] = ds[var].astype('int32')
    return ds

def save_to_netcdf(ds, output_file):
    """
    Save the xarray Dataset to a NetCDF file.
    """
    # Convert int64 to int32 for THREDDS compatibility
    ds = convert_int64_to_int32(ds)
    for attr in ['calendar', 'units']:
        if attr in ds['Timestamp'].attrs:
            del ds['Timestamp'].attrs[attr]

    ds.to_netcdf(output_file, mode='w', format='NETCDF4', engine='h5netcdf')

# %% ---- Process Data ---- #
glob_file_2025 =  DATA_PATH / "GLOB_data" / "GLOB_data_30sec_2025_NYA.dat"
bsrn_file_2025 =  DATA_PATH / "NYA_BSRN_data" / "NYA_radiation_2025-all.tab"
latitude=78.92240; longitude=11.92174 # Ny-Ålesund
f =10 # minute

# Load data
df = read_and_preprocess_data(glob_file_2025)
df_bsrn = pd.read_csv(bsrn_file_2025, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')

# Resample to f-minute intervals 
if f > 1: 
    freq = f'{f} minutes'
    df.index = df.index + pd.Timedelta(minutes=f/2)
    df = df.resample(f'{f}min').mean()
    df_bsrn = df_bsrn.resample(f'{f}min').first()
    output_file =  DATA_PATH / "GLOB_data" / f"GLOB_data_{f}min_2025.nc"
else: 
    freq = '30 seconds'
    output_file =  DATA_PATH / "GLOB_data" / "GLOB_data_30sec_2025.nc"

df_bsrn.index = df_bsrn.index.tz_localize('UTC')

df = substitute_ghi_ground(df, df_bsrn) # Substitute the ground reflection with BSRN data between 2025-05-09 and 2025-08-21

ds = convert_df_to_xarray_with_metadata(df, latitude, longitude)
ds = add_solar_angles_and_coordinates(ds, fct, latitude, longitude)
ds = compute_and_filter_albedo(ds)

ds = add_global_attributes(
ds=ds,
title= 'Global tilted irradiance on 25 planes measured in Ny-Ålesund (2025)',
summary=f'This file includes a time series of global tilted irradiance (GTI) \
measurements every {freq}, in 2025 in Ny-Ålesund.\n\
The measurements were acquired with a 25-pyranometer array called GLOB, mounted \
in the shape of rhombicuboctahedron. They consist of solar irradiance recorded \
with Apogee SP-110 pyranometers on 25 different planes.\n\
Each GTI components is labeled with the format "azimuth_tilt". We also included \
the solar angles associated to each timestamp, calculated with pvlib.',
keywords='GCMDSK: EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > SHORTWAVE RADIATION, '
'GCMDLOC: GEOGRAPHIC REGION > POLAR, \
GCMDLOC: CONTINENT > EUROPE > NORTHERN EUROPE > SCANDINAVIA > NORWAY',
geospatial_lat_min='78.92240',
geospatial_lat_max='78.92240',
geospatial_lon_min='11.92174',
geospatial_lon_max='11.92174',
time_coverage_start='2025-03-15T10:15:30',
time_coverage_end='2025-09-05T23:59:30',
history='- We created the file using netCDF4 in Python.',
comments='The albedo component was calculated with two upward and downward Apogee \
pyranometers from GLOB. But, between 2025-05-09 and 2025-08-22 the ground GHI \
sensor fell, and we used the albedo measurement from the BSRN station nearby instead.',
date_created= '2025-09-18',
creator_type='person, person, person',
creator_institution='The University Centre in Svalbard, \
The University Centre in Svalbard, The University Centre in Svalbard',
creator_name='Arthur Garreau, Aleksey Shestov, Sebastian Sikora',
creator_url='https://orcid.org/0000-0001-9509-1061, \
https://orcid.org/0000-0001-9601-8958, https://orcid.org/0009-0004-1874-7126',
creator_email='arthurg@unis.no, alekseys@unis.no, guliksen@gmail.com',
institution='The University Centre in Svalbard',
license='https://creativecommons.org/licenses/by/4.0/ (CC-BY-4.0)',
iso_topic_category='climatologyMeteorologyAtmosphere',
publisher_name='Arctic Data Centre',
publisher_institution='Norwegian Meteorological Institute',
publisher_url='https://adc.met.no',
publisher_email='adc-support@met.no',
station_name='GLOB Ny-Ålesund',
instrument_type='Apogee SP-110'
)
   
save_to_netcdf(ds, output_file)
print(f"{f}-minute NetCDF file created at: {output_file}")
# Print the dataset period
print('The ncdf dataset is from', str(ds['Timestamp'].values[0])[0:19], 'to'
      , str(ds['Timestamp'].values[-1])[0:19])

ds.close();