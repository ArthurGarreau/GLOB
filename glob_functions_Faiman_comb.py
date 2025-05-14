# -*- coding: utf-8 -*-
"""
Solar Angle and Irradiance Calculation Functions
================================================

This script provides a collection of functions for calculating solar angles, irradiance components,
and geometry coefficients using the Faiman et al. (1992) method. It includes functionalities for
estimating beam and diffuse irradiance, creating variable tables, and solving for irradiance components.

Key Features:
-------------
- Calculates solar angles (zenith, azimuth, declination, etc.) using pvlib.
- Computes the angle of incidence and zenith angle for a given plane.
- Estimates beam and diffuse irradiance components using the Faiman method.
- Provides geometry coefficients for beam and diffuse irradiance.
- Utilizes parallel processing for efficient computation of irradiance estimations.

Dependencies:
-------------
- numpy
- pvlib
- pandas
- joblib

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 6, 2024
"""

import numpy as np
import pvlib
import pandas as pd
from joblib import Parallel, delayed


def calculate_solar_angles(timestamps, latitude, longitude, altitude=6, temperature=-6):
    """
    Calculate solar angles (zenith, azimuth, declination, etc.) using pvlib.

    Parameters:
        timestamps (pd.DatetimeIndex): Timestamps for calculation.
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.
        altitude (float): Altitude of the location (default: 6 m).
        temperature (float): Temperature (default: -6°C).

    Returns:
        pd.DataFrame: Solar position data (zenith, azimuth, declination, etc.).
    """
    solar_position = pvlib.solarposition.get_solarposition(
        time=timestamps, latitude=latitude, longitude=longitude, altitude=altitude, temperature=temperature)
    
    eot = solar_position['equation_of_time']
    solar_position['hour_angle'] = pvlib.solarposition.hour_angle(timestamps, longitude, eot)
    solar_position['declination'] = np.degrees(pvlib.solarposition.declination_spencer71(timestamps.day_of_year))

    return solar_position

def incident_angle(solar_angles, plane_inclination, plane_azimuth, latitude, longitude):
    """
    Calculate the angle of incidence and zenith angle for a given plane.

    Parameters:
        solar_angles (pd.DataFrame): Dataframe of solar angles obtained with the above function "calculate_solar_angles" (degrees).
        plane_inclination (float): Inclination angle of the plane (degrees).
        plane_azimuth (float): Azimuth angle of the plane (degrees).
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.

    Returns:
        tuple: Angle of incidence and zenith angle (degrees).
    """
    
    delta = np.radians(solar_angles['declination'].values[0]) 
    omega = np.radians(solar_angles['hour_angle'].values[0])
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
    theta_i = np.degrees(theta_i)
    
    # zenith = solar_angles['zenith'].values[0]
    # azimuth = solar_angles['azimuth'].values[0]
    # theta_i = pvlib.irradiance.aoi(plane_inclination, plane_azimuth, zenith, azimuth) # in degree

    return theta_i

def Coef_b(solar_angles, plane_inclination, plane_azimuth, latitude, longitude, albedo):
    """
    Calculate the geometry coefficient for beam irradiance (Coef_b).

    Parameters:
        solar_angles (pd.DataFrame): Dataframe of solar angles obtained with the above function "calculate_solar_angles" (degrees).        plane_inclination (float): Inclination angle of the plane (degrees).
        plane_azimuth (float): Azimuth angle of the plane (degrees).
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.
        albedo (float): Surface albedo.

    Returns:
        float: Geometry coefficient for beam irradiance.
    """
    beta = plane_inclination #in degrees
    gamma = plane_azimuth #in degrees
    theta_i = incident_angle(solar_angles, beta, gamma, latitude, longitude)
    theta_z = solar_angles['zenith'].values[0]

    theta_i, theta_z, beta = np.radians(theta_i), np.radians(theta_z), np.radians(beta)
    
    # Calculate the sum of incidence and zenith angles
    Coef_b_value = np.cos(theta_i) + albedo * np.cos(theta_z)* 0.5 * (1 - np.cos(beta))
    
    return Coef_b_value


def Coef_d(plane_inclination, albedo):
    """
    Calculate the geometry coefficient for diffuse irradiance (Coef_d).

    Parameters:
        plane_inclination (float): Inclination angle of the plane (degrees).
        albedo (float): Surface albedo.

    Returns:
        float: Geometry coefficient for diffuse irradiance.
    """
    beta = np.radians(plane_inclination)
    Coef_d_value = 0.5 * (1 + np.cos(beta))*(1+np.sin(beta/2)**3) + albedo * 0.5 * (1 - np.cos(beta))
    return Coef_d_value


def estimation_diffuse_beam(variables, glob_value, solar_angles, lat, lon):
    """
    Estimate diffuse and beam irradiance using the Faiman method.

    Parameters:
        variables (list): List of GLOB plane names.
        glob_value (pd.DataFrame): Global irradiance values.
        solar_angles (pd.DataFrame): Dataframe of solar angles obtained with the above function "calculate_solar_angles" (degrees).
        lat (float): Latitude of the location.
        lon (float): Longitude of the location.

    Returns:
        tuple: Estimated diffuse and beam irradiance, and error.
    """
    timestamp = glob_value.index
    albedo = glob_value['albedo'].values[0]
    
    sigma = .05 #standard deviation error of pyrano
    
    # Generate all possible combinations of 2 variables
    table_azim_incli = create_variable_table(variables)
    vars_glob = table_azim_incli.index

    R_j = glob_value[vars_glob].values
    b_j = np.array( [Coef_b(solar_angles, 
                       table_azim_incli['inclination'].iloc[i], 
                       table_azim_incli['azimuth'].iloc[i], lat, lon, albedo) 
                   for i in range(len(table_azim_incli))] )
    d_j = np.array( [Coef_d(table_azim_incli['inclination'].iloc[i], albedo) 
                    for i in range(len(table_azim_incli))] )
    

    S_bb = np.sum(b_j**2/(sigma**2))
    S_dd = np.sum(d_j**2/(sigma**2))
    S_bd = np.sum(b_j*d_j/(sigma**2))
    S_Rd = np.sum(R_j*d_j/(sigma**2))
    S_Rb = np.sum(R_j*b_j/(sigma**2))
    delta = S_dd*S_bb - S_bd**2
    
    A_inter = np.array([[S_bb, -S_bd],
                        [-S_bd, S_dd]])
    A_moins_un = 1/delta * A_inter
    Y = np.array([S_Rd, S_Rb])
    # Estimate direct and diffuse components
    X_estim = np.matmul(A_moins_un, Y)
    D_prime, I_prime = X_estim
    
    # Calculate the variances = error
    D_prime_var = S_dd / delta ; I_prime_var = S_bb / delta

    # Calculate the solar position
    zenith = solar_angles['zenith'].values[0]
    cos_z = np.cos(np.radians(zenith))
    I_0 = pvlib.irradiance.get_extra_radiation(timestamp).values[0]
   
    # Solve D and I from D' and I' (Eq. 1.11 in Faiman et al. (1986))
    if not (np.isnan(X_estim[0]) and np.isnan(X_estim[1])):   
        D, I = solve_for_D_and_I(I_0, cos_z, D_prime, I_prime)
    
    # Reference values from ERBS model
    R_ghi = glob_value['GHI']
    ERBS = pvlib.irradiance.erbs(R_ghi, zenith, timestamp)
    I_ref = ERBS['dni']
    D_ref = ERBS['dhi']
        
    # Calculate the error in comparison to Erbs model
    error = ((I_ref - I)**2 + (D_ref - D)**2)**0.5
    # error = error.values[0]
    error = -delta

    if D > 0 and I > 0 and D_prime > 0 and I_prime > 0 and zenith < 90:
        return D, I, D_prime, I_prime, D_prime_var, I_prime_var, error
    else:
        # print("No valid combinations found (D and I might be negative or zenith<90°).")
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


# Optimized function to find the best combination
def find_best_combination(combs, glob_value, solar_angles, lat, lon):
    """
    Find the best combination of variables for irradiance estimation using a parralelisation method.

    Parameters:
        combs (list): List of variable combinations.
        glob_value (pd.DataFrame): Global irradiance values.
        zenith_angle (float): Solar zenith angle (degrees).
        lat (float): Latitude of the location.
        lon (float): Longitude of the location.

    Returns:
        tuple: Best estimated diffuse and beam irradiance, and the best combination.
    """
    # Evaluate all combinations in parallel
    results = Parallel(n_jobs=-1)(delayed(estimation_diffuse_beam)(
        list(comb), glob_value, solar_angles, lat, lon
    ) for comb in combs)

    # Extract errors from valid results
    errors = [result[6] for result in results]
    
    # Find the index of the minimum error among valid results
    best_index = np.nanargmin(errors)
    
    # Get the best result
    D, I, D_prime, I_prime, D_prime_var, I_prime_var, error = results[best_index]
    comb_opt = combs[best_index]

    return round(D), round(I), round(D_prime), round(I_prime), round(D_prime_var), round(I_prime_var), comb_opt

def create_variable_table(variables):
    """
    Create a table mapping GLOB planes to their azimuth and inclination.

    Parameters:
        variables (list): List of GLOB plane names.

    Returns:
        pd.DataFrame: Table with azimuth and inclination for each plane.
    """

    azimuth_mapping = {
        'S': 0, 'SW': 45, 'W': 90, 'NW': 135,
        'N': 180, 'NE': -135, 'E': -90, 'SE': -45
    }


    table = pd.DataFrame(columns=["azimuth", "inclination"])
    for var in variables:
        if var == 'GHI':
            table.loc[var] = [0, 0]
        else:
            direction, beta = var.split('_')
            table.loc[var] = [azimuth_mapping[direction], int(beta)]
    return table


def solve_for_D_and_I(I_0, cos_z, D_prime, I_prime):
    """
    Solve for D and I using the quadratic equation.

    Parameters:
        I_0 (float): Extraterrestrial irradiance.
        cos_z (float): Cosine of the zenith angle.
        D_prime (float): Intermediate value for diffuse irradiance.
        I_prime (float): Intermediate value for beam irradiance.

    Returns:
        tuple: Estimated diffuse and beam irradiance.
    """
    a = 1
    b = (I_0 - I_prime) * cos_z - D_prime
    c = -D_prime * I_0 * cos_z

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c

    # Check if the discriminant is non-negative
    if discriminant < 0:
        # print("The discriminant is negative, no real roots exist.")
        D, I = np.nan, np.nan
    else:
        # Calculate the two solutions for D
        D1 = (-b + np.sqrt(discriminant)) / (2*a)
        D2 = (-b - np.sqrt(discriminant)) / (2*a)
    
        # Select the positive root
        if D1 > 0:
            D = D1
            I = I_0 * (1 - D_prime / D)
        elif D2 > 0:
            D = D2
            I = I_0 * (1 - D_prime / D)
        else:
            # print("No positive root found for D.")
            D, I = np.nan, np.nan
          
    return D, I

def calc_Dprime_Iprime(I_0, zenith, D, I):
    """
    Calculate D' and I' from D and I using the provided formulas.

    Parameters:
        I_0 (float): Extraterrestrial irradiance.
        zenith (float): Zenith angle (in degree).
        D (float): Diffuse irradiance.
        I (float): Beam irradiance.

    Returns:
        tuple: Calculated D' and I'.
    """
    
    cos_z = np.cos(np.radians(zenith))
        
    # Calculate I' using the formula
    I_prime = I*(1 + D / (I_0 * cos_z))

    # Calculate D' using the formula
    D_prime = (I_0 - I) * D / I_0
    return D_prime, I_prime

