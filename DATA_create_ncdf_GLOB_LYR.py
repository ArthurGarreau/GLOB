# -*- coding: utf-8 -*-
"""
GLOB Data Processing and NetCDF Conversion Script
=================================================

This script processes GLOB data from a CSV file, calculates solar angles, computes surface albedo,
and converts the data into a NetCDF file format.

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

def resample_data(df, frequency):
    """
    Resample data to a specified frequency and compute the mean.
    """
    df.index = df.index + pd.Timedelta(minutes=frequency / 2)
    resampled_df = df.resample(f'{frequency}min').mean()
    return resampled_df

def convert_to_xarray_with_metadata(resampled_df, latitude, longitude):
    """
    Convert a resampled DataFrame to an xarray Dataset and add metadata to the Timestamp variable.

    Parameters:
        - resampled_df: Pandas DataFrame that has been resampled.

    Returns:
        - xarray Dataset with added metadata for the Timestamp variable.
    """
    ds = xr.Dataset.from_dataframe(resampled_df)
    timestamps_ds = pd.to_datetime(ds['Timestamp']).tz_localize('UTC')
    timestamps_ds = timestamps_ds.astype('int64') / 1e9

    ds['Timestamp'] = timestamps_ds.values
    ds['Timestamp'].attrs.update({
        'calendar': 'gregorian',
        'long_name': 'UTC time',
        'units':  'seconds since 1970-01-01T00:00:00Z',
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
            'long_name': 'Ground global horizontal irradiance',
            'standard_name': 'solar_irradiance',
            'units': 'W m-2',
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

    timestamps_ds = pd.to_datetime(ds['Timestamp'], unit='ns').tz_localize('UTC')
    solar_angles = fct.calculate_solar_angles(timestamps_ds, latitude, longitude)
    solar_angles_not_standard = solar_angles[['declination','equation_of_time', 'hour_angle', 'apparent_elevation', 'apparent_zenith']]
    solar_angles_standard = solar_angles.drop(columns=['declination','equation_of_time', 'hour_angle', 'apparent_elevation', 'apparent_zenith'])

    # Add solar angles to the dataset
    for var, values in solar_angles_standard.items():
        ds[var] = (('Timestamp',), values.values)
        ds[var].attrs.update({
            'units': 'degrees',
            'long_name': f'Solar {var.replace("_", " ")}',
            'standard_name': 'solar_' + var + '_angle'
        })
        
    for var, values in solar_angles_not_standard.items():
        ds[var] = (('Timestamp',), values.values)
        ds[var].attrs.update({
            'units': 'degrees',
            'long_name': f'Solar {var.replace("_", " ")}',
        })

    return ds

def compute_and_filter_albedo(ds, ds_KZ):
    """
    Compute albedo from the irradiance of the K&Z in Adventdalen, filter valid values, and calculate daily mean albedo between 10:00 and 12:00.
    Data ds_KZ are from :
        Garreau, A., Shestov, A., Sikora, S., & Sjöblom, A. (2024). Aggregated SW and LW irradiance downwelling and upwelling in Adventdalen [NetCDF4-CF]. Arctic Data Centre. https://doi.org/10.21343/psy9-3e97

    Parameters:
        - ds: xarray Dataset containing the data.

    Returns:
        - xarray Dataset with added albedo data.
    """
    ds_KZ = ds_KZ.rename({'time':'Timestamp'})
    albedo = ds_KZ['SWup'] / ds_KZ['SWdown']
    albedo = albedo.drop_duplicates(dim='Timestamp')
    albedo = albedo.where((albedo >= 0) & (albedo <= 1))
    filtered_albedo = albedo.where((albedo.indexes['Timestamp'].hour >= 10) & (albedo.indexes['Timestamp'].hour < 12))
    daily_mean_albedo = filtered_albedo.resample(Timestamp='1D').mean(skipna=True)
    
    # The instrument was out of order in July 2024 and therefore the albedo is set to the value it had just before this event. 
    # The ground was free of snow during this period so the assumption is realistic.
    mask = (daily_mean_albedo['Timestamp'] >= pd.Timestamp('2024-07-08')) & (daily_mean_albedo['Timestamp'] <= pd.Timestamp('2024-08-19'))
    daily_mean_albedo = daily_mean_albedo.where(~mask, 0.115) # Assumed value for the missing data period.

    ds_timestamp = pd.to_datetime(ds.indexes['Timestamp'], unit='s')
    aux_var = ds['GHI']; aux_var['Timestamp'] = ds_timestamp
    new_albedo = xr.full_like(aux_var, np.nan)
    for day in daily_mean_albedo.Timestamp.values:
        daily_mean = daily_mean_albedo.sel(Timestamp=day).item()
        new_albedo = new_albedo.where(new_albedo.Timestamp.dt.floor('D') != day, daily_mean)

    ds['albedo'] = (('Timestamp',), new_albedo.values)

    ds['albedo'].attrs.update({
        'units': '1',
        'long_name': 'Surface Albedo',
        'standard_name': 'surface_albedo'
    })

    return ds


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
        'keywords_vocabulary': kwargs.get('keywords_vocabulary', ''),
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
        'reference': kwargs.get('reference', ''),
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

    ds.to_netcdf(output_file, mode='w', format='NETCDF4', engine='h5netcdf')
    print(f"{f}-minute NetCDF file created at: {output_file}")


# %% ---- Process Data for each year ---- #

file_KZ = DATA_PATH / "Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"
latitude = 78.200318; longitude = 15.840308 # Adventdalen
f = .5  # minute

ds_all=[]
# Read and preprocess data for each year
for year in [2023, 2024]:
    glob_file = DATA_PATH / "GLOB_data" / f"GLOB_data_30sec_{year}.dat"
    # Resample to f-minute intervals 
    if f >= 1: 
        freq = f'{f} minutes'
        output_file =  DATA_PATH / "GLOB_data" / f"GLOB_data_{f}min_2023-24.nc"
        df = read_and_preprocess_data(glob_file)
        df = resample_data(df, f)
        ds = convert_to_xarray_with_metadata(df, latitude, longitude)
        ds_all.append(ds)
    else: 
        output_file =  DATA_PATH / "GLOB_data" / f"GLOB_data_{int(f*60)}sec_2023-24.nc"
        freq = f'{int(f*60)} seconds'
        df = read_and_preprocess_data(glob_file)
        ds = convert_to_xarray_with_metadata(df, latitude, longitude)
        ds_all.append(ds)
    print(year, "read and preprocessed")

# Merge datasets for 2023 and 2024
merged_ds = xr.merge(ds_all)

# Apply functions to the merged dataset
ds_KZ = fct.read_netcdf(file_KZ)
merged_ds = add_solar_angles_and_coordinates(merged_ds, fct, latitude, longitude)
merged_ds = compute_and_filter_albedo(merged_ds, ds_KZ)

merged_ds = add_global_attributes(
ds=merged_ds,
title= 'Global tilted irradiance on 25 planes measured in Adventdalen (2023-24)',
summary=f'This file includes a time series of global tilted irradiance (GTI) \
measurements every {freq}, from 2023 to 2024 in the valley of Adventdalen, Svalbard.\n\
The measurements were acquired with a 25-pyranometer array called GLOB, mounted \
in the shape of rhombicuboctahedron. They consist of solar irradiance recorded \
with Apogee SP-110 pyranometers on 25 different planes.\n\
Each GTI components is labeled with the format "azimuth_tilt". We also included \
the solar angles associated to each timestamp, calculated with pvlib.',
keywords='GCMDSK: EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > SHORTWAVE RADIATION, \
GCMDLOC: GEOGRAPHIC REGION > POLAR, \
GCMDLOC: CONTINENT > EUROPE > NORTHERN EUROPE > SCANDINAVIA > NORWAY',
keywords_vocabulary='GCMDSK: GCMD Science Keywords, GCMDLOC: GCMD Location',
geospatial_lat_min='78.200318',
geospatial_lat_max='78.200318',
geospatial_lon_min='15.840308',
geospatial_lon_max='15.840308',
time_coverage_start='2023-02-23T12:00:00',
time_coverage_end='2024-10-14T08:50:00',
history='- We created the file using netCDF4 in Python.',
comments='The albedo component was retrieved from the Kipp an Zonen CNR1 instrument \
100m away from GLOB.',
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
station_name='GLOB Adventdalen',
instrument_type='Apogee SP-110',
reference='https://doi.org/123/abc (Developped in Python. Citation: ...[to come])'
)

save_to_netcdf(merged_ds, output_file)
# Print the dataset period
print('The ncdf dataset is from ', str(merged_ds['Timestamp'].values[0])[0:19], 'to'
      , str(ds['Timestamp'].values[-1])[0:19])
ds.close(); merged_ds.close()

