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
        'long_name': 'Latitude of the location',
        'standard_name': 'latitude'
    })

    ds['longitude'].attrs.update({
        'units': 'degrees',
        'long_name': 'Longitude of the location',
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
    filtered_albedo = filtered_albedo.sortby('Timestamp')
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
    
    timestamps_ds = ds['Timestamp']
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
        
            theta_i = fct.incident_angle(solar_angles, inclinations, azimuths, lat, lon)
            theta_i[theta_i > 90] = 90; ####  !!!! VERY IMPORTANT !!!!
            # Calculate absolute air mass (AM_a)
            AM_a = np.exp(-0.0001184 * h) * (np.cos(np.radians(theta_i)) + 0.5057 * (96.080 - theta_i)**(-1.634))**(-1)

            # Calculate f1
            f1 = 0.000263 * (AM_a)**3 - 0.00632 * (AM_a)**2 + 0.054 * (AM_a) + 0.932

            # Calculate f2
            f2 = -0.00000004504 * theta_i**3 - 0.00001357 * theta_i**2 + 0.0006074 * theta_i + 1
            
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

# %% ---- Process Data ---- #
glob_file_2025 = DATA_PATH / "GLOB_data_30sec_2025_NYA.dat"
bsrn_file_2025 = DATA_PATH.parent / "Irradiance" / "NYA" / "NYA_radiation_2025-all.tab"
f = int(float(input("Please enter the data frequency required in minutes: "))) #minute
latitude=78.92240; longitude=11.92174 # Ny-Ålesund

# Load glob data
df = read_and_preprocess_data(glob_file_2025)
# Load the BSRN data
df_bsrn = pd.read_csv(bsrn_file_2025, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')

# Resample to f-minute intervals 
if f > 1: # minute
    df.index = df.index + pd.Timedelta(minutes=f/2)
    df = df.resample(f'{f}min').mean()
    
    df_bsrn = df_bsrn.resample(f'{f}min').first()
    freq = f'{f} minutes'
else: 
    freq = '30 seconds'

df_bsrn.index = df_bsrn.index.tz_localize('UTC')

df = substitute_ghi_ground(df, df_bsrn) # Substitute the ground reflection with BSRN data between 2025-05-09 and 2025-08-21

ds = convert_df_to_xarray_with_metadata(df)
ds = add_solar_angles_and_coordinates(ds, fct, latitude, longitude)
ds = compute_and_filter_albedo(ds)
ds = spectral_correction(ds)

ds = add_global_attributes(
    ds=ds,
    title= 'Tilted irradiance measurements on 25 orientations with GLOB in Ny-Ålesund (2025)',
    summary='This file includes a time series of global tilted irradiance (GTI) measurements every {freq}, in 2025 in Ny-Ålesund '
'Each GTI components is labeled with the format "azimuth_tilt". '
'The measurements consist of solar irradiance recorded with Apogee SP-110 pyranometers. '
'The pyranometers are mounted in a shape of a rhombicuboctahedron called "GLOB". '
'We included the solar angles associated to each timestamp, calculated with pvlib.',
    keywords='GCMDSK: EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > SHORTWAVE RADIATION, '
'GCMDLOC: GEOGRAPHIC REGION > POLAR, GCMDLOC: CONTINENT > EUROPE > NORTHERN EUROPE > SCANDINAVIA > NORWAY',
    geospatial_lat_min='78.92240',
    geospatial_lat_max='78.92240',
    geospatial_lon_min='11.92174',
    geospatial_lon_max='11.92174',
    time_coverage_start='2025-03-15',
    time_coverage_end='2025-09-05',
    history='- We created the file using netCDF4 in Python.\n'
'- We applied the spectral filter for silicon cell pyranometer data detailed in Balthazar et al. (2015)',
    comments='The albedo component was calculated with two Apogee from GLOB up and down. '
'But between 2025-05-09 and 2025-08-22 the ground GHI sensor fell, therefor we took albedo measurment from the BSRN station nearby using Kipp and Zonen instruments.',
    date_created= '2025-09-18',
    creator_type='person, person, person',
    creator_institution='The University Centre in Svalbard, The University Centre in Svalbard, The University Centre in Svalbard',
    creator_name='Arthur Garreau, Aleksey Shestov, Sebastian Sikora',
    creator_url='https://orcid.org/0000-0001-9509-1061, https://orcid.org/0000-0001-9601-8958, https://orcid.org/0009-0004-1874-7126',
    creator_email='arthurg@unis.no, alekseys@unis.no, guliksen@gmail.com',
    institution='The University Centre in Svalbard',
    license='https://creativecommons.org/licenses/by/4.0/ (CC-BY-4.0)',
    iso_topic_category='climatologyMeteorologyAtmosphere',
    publisher_name='Arctic Data Centre',
    publisher_institution='Norwegian Meteorological Institute',
    publisher_url='https://adc.met.no',
    publisher_email='adc-support@met.no',
    station_name='GLOB Ny-Ålesund',
    instrument_type='Apogee Silicon Cell SP-110',
)

if f>1: 
    output_file = DATA_PATH / f"GLOB_data_{f}min_2025.nc"
else: 
    output_file = DATA_PATH / "GLOB_data_30sec_2025.nc"
save_to_netcdf(ds, output_file)

# Print variables in the dataset
print('The ncdf dataset is from ', str(ds['Timestamp'].values[0])[0:10], 'to'
      , str(ds['Timestamp'].values[-1])[0:10])

