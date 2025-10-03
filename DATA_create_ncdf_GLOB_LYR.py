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
from config_path import DATA_PATH, GLOB_DATA_PATH
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
        'long_name': 'Latitude of the location',
        'standard_name': 'latitude'
    })

    ds['longitude'].attrs.update({
        'units': 'degrees',
        'long_name': 'Longitude of the location',
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
    
    timestamps_ds = ds['Timestamp']; timestamps_ds = pd.DatetimeIndex(timestamps_ds)
    lat = float(ds.latitude.values); lon = float(ds.longitude.values)

    # List of components to correct
    components = ['GHI', 'S_45', 'S_90', 'S_135', 'SW_45', 'SW_90', 'SW_135',
                  'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135', 'N_45',
                  'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45', 'E_90',
                  'E_135', 'SE_45', 'SE_90', 'SE_135']

    # Apply the correction to each component
    for component in components:
        if component in ds:

            table_azim_incli = fct.create_variable_table([component])
            inclinations = table_azim_incli['inclination'].values
            azimuths = table_azim_incli['azimuth'].values

            solar_angles = pd.DataFrame({
    'declination': ds['declination'].values,
    'hour_angle': ds['hour_angle'].values
})
        
            theta_i = fct.calculate_incident_angle(solar_angles, inclinations, azimuths, lat, lon)
            theta_i[theta_i > 90] = 90; ####  !!!! VERY IMPORTANT !!!!
            # Calculate absolute air mass (AM_a)
            AM_a = np.exp(-0.0001184 * h) * (np.cos(np.radians(theta_i)) + 0.5057 * (96.080 - theta_i)**(-1.634))**(-1)

            # Calculate f1
            f1 = 0.000263 * (AM_a)**3 - 0.00632 * (AM_a)**2 + 0.054 * (AM_a) + 0.932

            # Calculate f2
            f2 = -0.00000004504 * theta_i**3 - 0.00001357 * theta_i**2 + 0.0006074 * theta_i + 1
            
            ds[f'{component}'] = ds[component] / (f1 * f2)

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
        'keywords_vocabulary': kwargs.get('keywords_vocabulary', 'GCMDSK: GCMD Science Keywords, GCMDLOC: GCMD'),
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
    })

    # Add station_name and instrument_type as variables if provided
    if 'station_name' in kwargs:
        ds.coords['station_name'] = ((), kwargs['station_name'])
        ds['station_name'].attrs = {
            'long_name': "GLOB 25 pyranometers"
            }
    if 'instrument_type' in kwargs:
        ds.coords['instrument_type'] = ((), kwargs['instrument_type'])
        ds['instrument_type'].attrs = {
            'long_name': "25 pyranometers Apogee",
            } 

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

file_KZ = DATA_PATH.parent / "Irradiance_ncdf" / \
"Adventdalen_global_horizontal_irradiances_LW_SW_all.nc"
latitude = 78.200318; longitude = 15.840308 # Adventdalen
f = .5  # minute

ds_all=[]
# Read and preprocess data for each year
for year in [2023, 2024]:
    # Resample to f-minute intervals 
    if f >= 1: 
        file_path = GLOB_DATA_PATH / f"GLOB_data_{f}min_{year}.dat"; freq = f'{f} minutes'
        output_file = GLOB_DATA_PATH / f"GLOB_data_{f}min_2023-24.nc"
        df = read_and_preprocess_data(file_path)
        df.index = df.index + pd.Timedelta(minutes=f/2)
        df = df.resample(f'{f}min').mean()
        ds = convert_to_xarray_with_metadata(df)
        ds_all.append(ds)
    else: 
        file_path = GLOB_DATA_PATH / f"GLOB_data_{int(f*60)}sec_{year}.dat"
        output_file = GLOB_DATA_PATH / f"GLOB_data_{int(f*60)}sec_2023-24.nc"
        freq = f'{int(f*60)} seconds'
        df = read_and_preprocess_data(file_path)
        ds = convert_to_xarray_with_metadata(df)
        ds_all.append(ds)
    print(year, "read and preprocessed")

# Merge datasets for 2023 and 2024
merged_ds = xr.merge(ds_all)

# Apply functions to the merged dataset
ds_kZ = xr.open_dataset(file_KZ)
merged_ds = add_solar_angles_and_coordinates(merged_ds, fct, latitude, longitude)
merged_ds = compute_and_filter_albedo(merged_ds, ds_kZ)
merged_ds = spectral_correction(merged_ds)

merged_ds = add_global_attributes(
ds=merged_ds,
title= 'Global tilted irradiance on 25 faces measured in Adventdalen (2023-24)',
summary=f'This file includes a time series of global tilted irradiance (GTI) \
measurements every {freq}, from 2023 to 2024 in the valley of Adventdalen, Svalbard.\n\
The measurements were acquired with a 25-pyranometer array called GLOB, mounted \
in the shape of rhombicuboctahedron. They consist of solar irradiance recorded \
with Apogee SP-110 pyranometers.\n\
Each GTI components is labeled with the format "azimuth_tilt". We also included \
the solar angles associated to each timestamp, calculated with pvlib.',
keywords='GCMDSK: EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > SHORTWAVE RADIATION, '
'GCMDLOC: GEOGRAPHIC REGION > POLAR, \
GCMDLOC: CONTINENT > EUROPE > NORTHERN EUROPE > SCANDINAVIA > NORWAY',
geospatial_lat_min='78.200318',
geospatial_lat_max='78.200318',
geospatial_lon_min='15.840308',
geospatial_lon_max='15.840308',
time_coverage_start='2023-02-23',
time_coverage_end='2024-10-14',
history='- We created the file using netCDF4 in Python.\n'
'- We applied the spectral filter for silicon cell pyranometer data detailed in \
Balthazar et al. (2015)',
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
instrument_type='Apogee SP-110'
)

save_to_netcdf(merged_ds, output_file)
ds.close(); ds_kZ.close(); merged_ds.close()
