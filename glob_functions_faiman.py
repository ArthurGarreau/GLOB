# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 14:33:44 2024

@author: arthurg

Contact: arthurg@unis.no
"""


# %% FUNCTION FOR GLOB CALCULATIONS %%


import numpy as np
import pvlib
import pandas as pd
from itertools import combinations


def calculate_solar_angles(timestamps, latitude, longitude, altitude=6, temperature=-6):
    """
    Calculate solar angles using pvlib for the given timestamps and location.
    """

    solar_position = pvlib.solarposition.get_solarposition(
        time=timestamps, latitude=latitude, longitude=longitude, altitude=altitude, temperature=temperature
    )
   
    
    eot = solar_position['equation_of_time']
    solar_position['hour_angle'] = pvlib.solarposition.hour_angle(timestamps, longitude, eot)
    solar_position['declination'] = pvlib.solarposition.declination_spencer71(timestamps.day_of_year)
    
    return solar_position

def incident_and_zenith_angle(timestamps, plane_inclination, plane_azimuth, latitude, longitude, altitude=6, temperature=-6):
    """
    Parameters:
    - timestamps: pd.DatetimeIndex, the timestamps to calculate angles for.
    - angles: all in degree
    
    Calculate the cosine of the angle of incidence, θ.
    """
    solar_angles = calculate_solar_angles(timestamps, latitude, longitude)
    
    delta = np.radians(solar_angles['declination'])
    omega = np.radians(solar_angles['hour_angle'])
    beta = np.radians(plane_inclination)
    gamma = np.radians(plane_azimuth)
    phi = np.radians(latitude)
    
    theta_i = np.arccos(
        np.sin(delta) * np.sin(phi) * np.cos(beta) 
        - np.sin(delta) * np.cos(phi) * np.sin(beta) * np.cos(gamma) 
        + np.cos(delta) * np.cos(phi) * np.cos(beta) * np.cos(omega) 
        + np.cos(delta) * np.sin(phi) * np.sin(beta) * np.cos(omega) * np.cos(gamma) 
        + np.cos(delta) * np.sin(beta) * np.sin(omega) * np.sin(gamma)
    )
    
    theta_z = solar_angles['apparent_zenith'].values
    
    return np.degrees(theta_i), theta_z

def C_b(timestamps, plane_inclination, plane_azimuth, latitude, longitude, albedo, altitude=6, temperature=-6):
    """
    Calculate the sum of the incidence angle and the zenith angle for given timestamps and location.

    Parameters:
    - timestamps: pd.DatetimeIndex, the timestamps to calculate angles for.
    - latitude: float, latitude of the location in degrees.
    - longitude: float, longitude of the location in degrees.
    - beta: float, tilt angle of the surface in degrees.
    - gamma: float, azimuth angle of the surface in degrees.
    - altitude: float, altitude of the location in meters (default is 6 m).
    - temperature: float, ambient temperature in Celsius (default is -6 °C).

    Returns:
    - C_b_values: np.ndarray, the sum of incidence and zenith angles for each timestamp.
    """

    beta = plane_inclination
    gamma = plane_azimuth
    # Calculate incidence angle
    theta_i, theta_zenith = incident_and_zenith_angle(timestamps, beta, gamma, latitude, longitude)

    theta_i, theta_zenith = np.radians(theta_i), np.radians(theta_zenith)
    theta_zenith[theta_zenith < 0 ] = np.nan; theta_i[theta_i < 0 ] = np.nan; 
    
    beta = np.radians(plane_inclination)

    # Calculate the sum of incidence and zenith angles
    C_b_value = np.cos(theta_i) + albedo * np.cos(theta_zenith)* 0.5 * (1 - np.cos(beta))
    return C_b_value.iloc[0]


def C_d(plane_inclination, albedo):
    """
    Calculate the geometry coefficient for Diffuse irradiance (C_{d}).
    Assumes isotropic diffuse model.
    """
    beta = np.radians(plane_inclination)
    C_d_value = 0.5 * (1 + np.cos(beta)) + albedo * 0.5 * (1 - np.cos(beta))
    return C_d_value.iloc[0] 


def find_best_estimation(combinations_2, glob_value, NYA_value, orientations_dict, ds_glob):
    """
    Find the best combination of orientation indices that minimizes the error in estimating direct
    and diffuse irradiance components.

    Parameters:
        combinations_2 (list): List of tuples containing index pairs.
        glob_value (dict): Dictionary with global irradiance values.
        NYA_value (dict): Dictionary with true direct and diffuse irradiance values.
        orientations_dict (dict): Dictionary with azimuth and tilt (beta) values for orientations.
        ds_glob (xarray.Dataset): Dataset containing latitude and other solar data.

    Returns:
        tuple: Best combination of orientations and the corresponding estimated components.
    """
    min_error = float('inf')
    best_combination = None
    best_X_estim = None

    # True values of direct (DIR) and diffuse (DIF) irradiance
    X_truth = np.array([NYA_value['DIR [W/m**2]'], NYA_value['DIF [W/m**2]']])
    
    timestamps = glob_value.index
    albedo = glob_value['albedo']
    latitude = ds_glob.latitude.values; longitude = ds_glob.longitude.values
    
    for comb in combinations_2:
        i, j = comb
        
        
        
        # Measured global irradiance values for the selected orientations
        Y = np.array([glob_value[i], glob_value[j]])

        # Calculate coefficients for the selected orientations
        b_i =  C_b(timestamps, 
                   orientations_dict[i]['Beta'], orientations_dict[i]['Azimuth'], 
                   latitude, longitude, albedo)

        b_j =  C_b(timestamps, 
                   orientations_dict[j]['Beta'], orientations_dict[j]['Azimuth'], 
                   latitude, longitude, albedo)
        
        d_i = C_d(orientations_dict[i]['Beta'], albedo)
        d_j = C_d(orientations_dict[j]['Beta'], albedo)

        # Create the matrix A and its inverse
        A = np.array([[b_i, d_i],
                      [b_j, d_j]])
        A_moins_un = np.linalg.inv(A)
        
        # print(A, 'and A-1 = ', A_moins_un)
        
        # Estimate direct and diffuse components
        X_estim = np.matmul(A_moins_un, Y)

        # Calculate the error (mean absolute error)
        error = np.mean(np.sqrt((X_truth - X_estim)**2))

        # Update the best result if the current error is smaller
        if error < min_error:
            min_error = error
            best_combination = comb
            best_X_estim = X_estim
            A_opt = A_moins_un

    # Print the current estimation vs truth for analysis
    print(f"Timestamp: {timestamps[0]} Combination {best_combination}: Estimated = {np.squeeze(best_X_estim)}, Truth = {np.squeeze(X_truth)}, Error = {np.squeeze(min_error)}")

    return best_combination, best_X_estim

def create_variable_dict(variables):
    """
    Create a dictionary structure mapping variables to their azimuth and beta values.

    Parameters:
        variables (list): List of variable names.

    Returns:
        dict: A dictionary with variable names as keys and their attributes as values.
    """
    azimuth_mapping = {
        'S': 0,
        'SW': 45,
        'W': 90,
        'NW': 135,
        'N': 180,
        'NE': -45,
        'E': -90,
        'SE': -135
    }

    variable_dict = {}

    for var in variables:
        if var == 'GHI':
            variable_dict[var] = {'Azimuth': 0, 'Beta': 0}
            continue

        parts = var.split('_')
        direction = parts[0]
        beta = int(parts[1])
        azimuth = azimuth_mapping.get(direction, None)

        variable_dict[var] = {'Azimuth': azimuth, 'Beta': beta}

    return variable_dict