# -*- coding: utf-8 -*-
"""
Solar Angle and Irradiance Calculation Functions
================================================
This script provides a collection of functions for calculating solar angles, irradiance components,
and geometry coefficients. The beam and diffuse estimation are based on the method described in Faiman et al. (1987).
Faiman, D., Zemel, A., & Zangvil, A. (1987). A method for monitoring insolation
in remote regions. Solar Energy, 38(5), 327–333. https://doi.org/10.1016/0038-092X(87)90004-1

Key Features:
-------------
- Calculates solar angles (zenith, azimuth, declination, etc.) using pvlib.
- Computes the angle of incidence for a given plane.
- Estimates beam and diffuse irradiance components using the Faiman method.
- Provides geometry coefficients for beam, diffuse, and reflected irradiance.
- Utilizes parallel processing for efficient computation of irradiance estimations.

Dependencies:
-------------
- numpy
- pvlib
- pandas
- scipy
- joblib
- re

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 6, 2024
"""

import numpy as np
import pvlib
import pandas as pd
from scipy.optimize import least_squares
from joblib import Parallel, delayed
import re

def calculate_solar_angles(timestamps, latitude, longitude, altitude=6, temperature=-6):
    """
    Calculate solar angles (zenith, azimuth, declination, etc.) using pvlib.

    Parameters:
        timestamps (pd.DatetimeIndex or array-like): Timestamps for calculation.
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.
        altitude (float): Altitude of the location (default: 6 m).
        temperature (float): Temperature (default: -6°C).

    Returns:
        pd.DataFrame: Solar position data (zenith, azimuth, declination, etc.).
    """
    if isinstance(timestamps, pd.core.indexes.datetimes.DatetimeIndex | pd._libs.tslibs.timestamps.Timestamp):
        if timestamps.tzinfo is None:
            timestamps = timestamps.tz_localize('UTC')
        solar_position = pvlib.solarposition.get_solarposition(
            time=timestamps, latitude=latitude, longitude=longitude, altitude=altitude, temperature=temperature)
        eot = solar_position['equation_of_time']
        solar_position['hour_angle'] = pvlib.solarposition.hour_angle(
            timestamps, longitude, eot)  # timestamps is already a pd.DatetimeIndex
    else:
        timestamps = pd.DatetimeIndex(timestamps)
        if timestamps.tzinfo is None:
            timestamps = timestamps.tz_localize('UTC')
        solar_position = pvlib.solarposition.get_solarposition(
            time=timestamps, latitude=latitude, longitude=longitude, altitude=altitude, temperature=temperature)
        eot = solar_position['equation_of_time']
        solar_position['hour_angle'] = pvlib.solarposition.hour_angle(
            timestamps, longitude, eot)  # timestamps is not a pd.DatetimeIndex

    solar_position['declination'] = np.degrees(pvlib.solarposition.declination_spencer71(timestamps.day_of_year))
    return solar_position

def calculate_incident_angle(solar_angles, plane_inclination, plane_azimuth, latitude, longitude):
    """
    Calculate the angle of incidence of the sun for a given plane.

    Parameters:
        solar_angles (pd.DataFrame): Dataframe or object of solar angles 
        requiring 'declination' and 'hour_angle' (degrees).
        plane_inclination (float): Inclination angle of the plane (degrees).
        plane_azimuth (float): Azimuth angle of the plane (degrees).
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.

    Returns:
        np.ndarray: Angle of incidence (degrees).
    """
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
    )  # see formula in Duffie, J. A., & Beckman, W. A. (2013).
       # Available Solar Radiation Ch. 1. In Solar Engineering of Thermal Processes (pp. 43–137).
       # John Wiley & Sons, Inc. https://doi.org/10.1002/9781118671603.ch2
    theta_i = np.degrees(theta_i)
    return theta_i

def calculate_beam_coefficient(solar_angles, plane_inclination, plane_azimuth, latitude, longitude):
    """
    Calculate the transposition coefficient for beam irradiance.

    Parameters:
        solar_angles (pd.DataFrame): Dataframe or object of solar angles 
        requiring 'zenith', 'declination' and 'hour_angle' (degrees).
        plane_inclination (float): Inclination angle of the plane (degrees).
        plane_azimuth (float): Azimuth angle of the plane (degrees).
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.

    Returns:
        np.ndarray: Beam coefficient.
    """
    theta_i = calculate_incident_angle(solar_angles, plane_inclination, plane_azimuth, latitude, longitude)
    theta_i[theta_i > 90] = 90  # !!! VERY IMPORTANT !!! The incident angles above 90 have to be set to 90 to avoid the beam having an influence
    theta_i = np.radians(theta_i)
    theta_z = np.radians(solar_angles['zenith'])
    beam_coefficient = np.cos(theta_i) / np.cos(theta_z)
    return beam_coefficient

def calculate_diffuse_coefficient(plane_inclination):
    """
    Calculate the transposition coefficient for diffuse irradiance.

    Parameters:
        plane_inclination (float): Inclination angle of the plane (degrees).

    Returns:
        np.ndarray: Diffuse coefficient.
    """
    beta = np.radians(plane_inclination)
    diffuse_coefficient = 0.5 * (1 + np.cos(beta))
    return diffuse_coefficient

def calculate_reflected_coefficient(plane_inclination, plane_azimuth=None, solar_angles=None):
    """
    Calculate the transposition coefficient for reflected irradiance.

    Parameters:
        plane_inclination (float): Inclination angle of the plane (degrees).

    Returns:
        np.ndarray: Reflected coefficient.
    """
    beta = np.radians(plane_inclination)
    reflected_coefficient = 0.5 * (1 - np.cos(beta))
    return reflected_coefficient

def least_squares_residuals(params, b, d, r, GTI, albedo=None):
    """
    Calculate the residuals for the least squares estimation of the quantity:
        GTI - (b + r*Albedo) * Beam_prime + (d + r*Albedo) * Diffuse_prime

    Parameters:
        params (list or array-like): The parameters to estimate
            [Diffuse_prime, Beam_prime, Albedo] if albedo==None (this is the nonlinear case).
            [Diffuse_prime, Beam_prime] if albedo is input (this is the linear case).
        b, d, r (array-like): The transposition coefficients for Beam, Diffuse, and Reflected.
        GTI (array-like): The observed values of GTI.

    Returns:
        np.ndarray: The residuals.
    """
    if albedo is None:  # nonlinear method case
        x, y, z = params  # (x,y,z) = (diffuse, beam, albedo)
        GTI_estim = (b + r * z) * y + (d + r * z) * x
    else:  # linear method case
        x, y = params  # (x,y) = (diffuse, beam)
        GTI_estim = (b + r * albedo) * y + (d + r * albedo) * x
    return GTI_estim - GTI

def estimate_diffuse_beam_faiman(variables, glob_value, lat, lon, method='linear'):
    """
    Estimate diffuse and beam irradiance based on either linear or nonlinear equations.

    Parameters:
        variables (list): List of GLOB plane names.
        glob_value (pd.DataFrame): Global irradiance and solar angles values
        obtained from the NetCDF GLOB dataset.
        lat (float): Latitude of the location.
        lon (float): Longitude of the location.
        method (str): Method to use for estimation ('linear' or 'nonlinear').

    Returns:
        list: Estimated diffuse and beam irradiance, albedo, and error.
    """
    
    table_azim_incli = create_variable_table(variables)
    GTI_glob = np.asarray(glob_value[table_azim_incli.index], dtype=float)

    # Interpolate NaN values in GTI_glob
    inclinations = table_azim_incli['inclination'].values
    azimuths = table_azim_incli['azimuth'].values
    
    # Necessary solar angles for the geometry calculation.
    # Can also be calculated with the function 'solar_angle_calculation'
    # if no in the initial dataset.
    solar_angles = glob_value[['zenith','hour_angle', 'declination']] 
    
    b = np.array(calculate_beam_coefficient(solar_angles, inclinations, azimuths, lat, lon))
    d = np.array(calculate_diffuse_coefficient(inclinations))
    r = np.array(calculate_reflected_coefficient(inclinations))
    X = np.column_stack((b, d, r))
    # Remove the NaN of the coefficients table
    rows_with_nan = np.isnan(X).any(axis=1)
    X = X[~rows_with_nan]
    b, d, r = X[:, 0], X[:, 1], X[:, 2]
    try:
        if method == 'linear':
            albedo = glob_value['albedo']
            initial_guess = [10, 10]
            # Add bounds if applicable
            bounds = ([0, 0], [np.inf, np.inf])
            results = least_squares(
                least_squares_residuals,
                initial_guess,
                bounds=bounds,
                args=(b, d, r, GTI_glob, albedo),
                max_nfev=1000  # Increase maximum number of function evaluations
            )
            D_prime, B_prime = results.x
        elif method == 'nonlinear':
            initial_guess = [10, 10, 10]
            # Add bounds if applicable
            bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
            results = least_squares(
                least_squares_residuals,
                initial_guess,
                bounds=bounds,
                args=(b, d, r, GTI_glob),
                max_nfev=1000  # Increase maximum number of function evaluations
            )
            D_prime, B_prime, albedo = results.x
        # Check if the solution is valid
        if not results.success:
            print("The least squares did not converge.")
            return [np.nan, np.nan, np.nan, np.nan, np.nan]
    except Exception:
        return [np.nan, np.nan, np.nan, np.nan, np.nan]
    
    zenith = glob_value['zenith']  # degrees
    cos_z = np.cos(np.radians(zenith))
    I_0 = pvlib.irradiance.get_extra_radiation(glob_value.name) # glob_value.name = timestamp
    if not (np.isnan(D_prime) and np.isnan(B_prime)):
        D, B = calculate_D_and_B(I_0, cos_z, D_prime, B_prime)
    D = D[0] if np.shape(D) == (1,) and not np.isnan(D[0]) else D
    B = B[0] if np.shape(B) == (1,) and not np.isnan(B[0]) else B
    D_prime = D_prime[0] if np.shape(D_prime) == (1,) and not np.isnan(D_prime[0]) else D_prime
    B_prime = B_prime[0] if np.shape(B_prime) == (1,) and not np.isnan(B_prime[0]) else B_prime
    albedo = albedo[0] if np.shape(albedo) == (1,) and not np.isnan(albedo[0]) else albedo
    if D > 0 and B > 0 and zenith < 90:
        return [np.round(D), np.round(B), np.round(albedo, 2), np.round(D_prime), np.round(B_prime)]
    else:
        return [np.nan, np.nan, np.nan, np.nan, np.nan]

def estimate_diffuse_beam_monte_carlo(variables, glob_value, lat, lon, n_simulations, error):
    """
    Modified estimate_diffuse_beam_faiman function for performing the Monte Carlo calculations.
    This function is made for calculating the error propagation of the 
    measurements in the estimations.

    Parameters:
        variables (list): List of GLOB plane names.
        glob_value (pd.DataFrame): Global irradiance and solar angles values
        obtained from the NetCDF GLOB dataset.
        lat (float): Latitude of the location.
        lon (float): Longitude of the location.
        n_simulations (int): Number of Monte Carlo simulations.
        error (float): Error margin for simulations.

    Returns:
        tuple: Estimated diffuse and beam irradiance.
    """
    table_azim_incli = create_variable_table(variables)
    GTI_glob = np.asarray(glob_value[table_azim_incli.index], dtype=float)
   
    inclinations = table_azim_incli['inclination'].values
    azimuths = table_azim_incli['azimuth'].values
    
    # Necessary solar angles for the geometry calculation.
    # Can also be calculated with the function 'solar_angle_calculation'
    # if no in the initial dataset.
    solar_angles = glob_value[['zenith','hour_angle', 'declination']] # necessary solar angles for the geometry calculation.

    b = np.array(calculate_beam_coefficient(solar_angles, inclinations, azimuths, lat, lon))
    d = np.array(calculate_diffuse_coefficient(inclinations))
    r = np.array(calculate_reflected_coefficient(inclinations))
    X = np.column_stack((b, d, r))
    # Remove the NaN of the coefficients table
    rows_with_nan = np.isnan(X).any(axis=1)
    X = X[~rows_with_nan]
    b, d, r = X[:, 0], X[:, 1], X[:, 2]

    try:
        # Monte Carlo simulation
        B_estimates = np.zeros(n_simulations)
        D_estimates = np.zeros(n_simulations)
        albedo = glob_value['albedo']

        for i in range(n_simulations):
            albedo_noisy = albedo * np.random.normal(1, error, size=1)
            initial_guess = [1, 1]
            # Add bounds if applicable
            bounds = ([0, 0], [np.inf, np.inf])
            
            results = least_squares(
                least_squares_residuals,
                initial_guess,
                bounds=bounds,
                args=(b, d, r, GTI_glob, albedo_noisy),
                max_nfev=1000  # Increase maximum number of function evaluations
            )
            D_prime, B_prime = results.x
            zenith = glob_value['zenith']  # degrees
            cos_z = np.cos(np.radians(zenith))
            I_0 = pvlib.irradiance.get_extra_radiation(glob_value.name) # glob_value.name = timestamp

            if not (np.isnan(D_prime) and np.isnan(B_prime)):
                D, B = calculate_D_and_B(I_0, cos_z, D_prime, B_prime)
                if D > 1 and D < 1372 and B > 1 and B < 1372:
                    D_estimates[i] = D
                    B_estimates[i] = B        
                else:
                    D_estimates[i] = np.nan
                    B_estimates[i] = np.nan
                    
    except Exception:
        return (np.nan, np.nan)
    if zenith < 90:
        return (D_estimates, B_estimates)
    else:
        return (np.nan, np.nan)

def find_best_estimation(combs, glob_value, lat, lon, true_estimation, method='linear'):
    """
    Find the best combination of variables for the least-square estimation of
    beam and diffuse (and albedo) using a parallelization method.

    Parameters:
        combs (list): List of variable combinations.
        glob_value (pd.DataFrame): Global irradiance and solar angles values
        obtained from the NetCDF GLOB dataset.
        lat (float): Latitude of the location.
        lon (float): Longitude of the location.
        true_estimation (list of lists): List of measured [diffuse, beam].

    Returns:
        list: Best estimated diffuse and beam irradiance, albedo, and the best combination.
    """
    # Evaluate all combinations in parallel
    results = np.array(
        Parallel(n_jobs=-1)(delayed(estimate_diffuse_beam_faiman)(
            list(comb), glob_value, lat, lon, method
        ) for comb in combs)
    )

    true_diffuse, true_beam = true_estimation[0], true_estimation[1]
    # Extract errors from valid results
    errors = [np.sqrt((results[:, 0] - true_diffuse)**2 + (results[:, 1] - true_beam)**2)]

    # Check if all results are NaN
    if np.all(np.isnan(errors)):
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]  # We add a Nan for the pyrano combination

    else:
        # Find the index of the minimum error among valid results
        best_index = np.nanargmin(errors)

        # Get the best result
        D, B, albedo, D_prime, B_prime = results[best_index, :]
        comb_opt = combs[best_index]

        return [np.round(D), np.round(B), np.round(albedo, 2), np.round(D_prime), np.round(B_prime), comb_opt]


def calculate_D_and_B(I_0, cos_z, D_prime, B_prime):
    """
    Solve for D and B using the quadratic equation.
    See equation 1.11 (Appendix A) in Faiman et al., (1987).

    Parameters:
        I_0 (float): Extraterrestrial irradiance.
        cos_z (float): Cosine of the zenith angle.
        D_prime (float): Intermediate value for diffuse irradiance.
        B_prime (float): Intermediate value for beam irradiance.

    Returns:
        tuple: Estimated diffuse and beam irradiance.
    """
    a = 1
    b = (I_0 - B_prime) * cos_z - D_prime
    c = -D_prime * I_0 * cos_z
    # Calculate the discriminant
    discriminant = b**2 - 4*a*c
    # Check if the discriminant is non-negative
    if discriminant < 0:
        D, B = np.nan, np.nan
    else:
        # Calculate the two solutions for D
        D1 = (-b + np.sqrt(discriminant)) / (2*a)
        D2 = (-b - np.sqrt(discriminant)) / (2*a)
        # Select the positive root
        if D1 > 0:
            D = D1
            B = I_0 / cos_z * (1 - D_prime / D)
        elif D2 > 0:
            D = D2
            B = I_0 / cos_z * (1 - D_prime / D)
        else:
            D, B = np.nan, np.nan
    return D, B

def calculate_Dprime_and_Bprime(I_0, zenith, D, I):
    """
    Calculate D' and I' from D and I using the provided formulas.
    See equation 1.9 (Appendix A) in Faiman et al., (1987).

    Parameters:
        I_0 (float): Extraterrestrial irradiance.
        zenith (float): Zenith angle (in degrees).
        D (float): Diffuse irradiance.
        I (float): Beam irradiance.

    Returns:
        tuple: Calculated D' and I'.
    """
    cos_z = np.cos(np.radians(zenith))
    # Calculate I' using the formula
    I_prime = I * (cos_z + D / I_0)
    # Calculate D' using the formula
    D_prime = (I_0 - I * cos_z) * D / I_0
    return D_prime, I_prime

def calculate_GTI_for_orientations(
    solar_angles, tilt_angles, azimuth_angles_calc, azimuth_angles_names,
    beam_prime, diffuse_prime, albedo, theta_z, lat, lon):
    """
    Calculate Global Tilted Irradiance for each combination of tilt and azimuth angles.

    Parameters:
    -----------
    solar_angles : DataFrame
        Solar angles containing zenith, azimuth, hour_angle, and delcination.
    tilt_angles : array-like
        Array of tilt angles.
    azimuth_angles_calc : array-like
        Array of azimuth angles for calculation.
    azimuth_angles_names : array-like
        Array of azimuth angles for column naming.
    beam_prime : array-like
        Beam irradiance values.
    diffuse_prime : array-like
        Diffuse irradiance values.
    albedo : array-like
        Albedo values.
    theta_z : array-like
        Zenith angles.
    lat : float
        Latitude of the location.
    lon : float
        Longitude of the location.

    Returns:
    -------
    df_results : DataFrame
        DataFrame with irradiance for each tilt and azimuth combination.
    """
    # Initialize a list to store results for each timestamp
    column_names = [f'gti{orientation}_{tilt}' for orientation in azimuth_angles_names for tilt in tilt_angles]
    df_results = pd.DataFrame(columns=column_names)

    for tilt in tilt_angles:
        for idx, azimuth in enumerate(azimuth_angles_calc):
            # Calculate incident angle
            theta_i = calculate_incident_angle(solar_angles, tilt, azimuth, lat, lon)
            theta_i[theta_i > 90] = 90

            # Calculate components and irradiance
            b = np.cos(np.radians(theta_i)) / np.cos(np.radians(theta_z))
            d = (1 + np.cos(np.radians(tilt))) / 2
            r = (1 - np.cos(np.radians(tilt))) / 2

            R = (albedo * r + b) * beam_prime + (albedo * r + d) * diffuse_prime

            # Store results
            df_results[f'gti{azimuth_angles_names[idx]}_{tilt}'] = R

    return df_results

def azimuth_to_orientation(azimuth):
    """
    Convert an azimuth angle to its corresponding cardinal orientation.

    Parameters:
        azimuth (float): The azimuth angle in degrees, which can be positive or negative.

    Returns:
        str: The cardinal orientation corresponding to the azimuth angle.
    """
    # Determine the orientation based on the azimuth
    if azimuth == 0:
        return 'N'
    elif azimuth == 45:
        return 'NE'
    elif azimuth == 90:
        return 'E'
    elif azimuth == 135:
        return 'SE'
    elif azimuth == 180 or azimuth == -180:
        return 'S'
    elif azimuth == 225:
        return 'SW'
    elif azimuth == 270:
        return 'W'
    elif azimuth == 315:
        return 'NW'
    else:
        return 'Unknown'

def orientation_to_azimuth(orientation):
    """
    Convert a cardinal orientation to its corresponding azimuth angle.

    Parameters:
        orientation (str): The cardinal orientation.

    Returns:
        float: The azimuth angle in degrees corresponding to the cardinal orientation.
    """
    # Define the azimuth mapping
    azimuth_mapping = {
        'S': 0,
        'SW': 45,
        'W': 90,
        'NW': 135,
        'N': 180,
        'NE': -135,
        'E': -90,
        'SE': -45
    }
    # Retrieve the azimuth angle based on the orientation
    return azimuth_mapping.get(orientation)

def find_missing_elements(input_list):
    """
    Find missing elements from the reference list of pyranometer names.

    Parameters:
        input_list (list): List of pyranometer names.

    Returns:
        list: Missing elements from the reference list.
    """
    reference_list = [
        'GHI', 'N_45', 'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45', 'E_90',
        'E_135', 'SE_45', 'SE_90', 'SE_135', 'S_45', 'S_90', 'S_135', 'SW_45',
        'SW_90', 'SW_135', 'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135'
    ]
    # Convert both lists to sets
    input_set = set(input_list)
    # Find the difference between the reference set and the input set
    missing_elements = [item for item in reference_list if item not in input_set]
    return missing_elements

def find_content_between_braces(file_path, line_number):
    """
    Function to return the name of the pyranometers that have not been used for
    the estimations of GTI.

    Parameters:
        file_path (str): Path to the file.
        line_number (int): Line number to search.

    Returns:
        list: Names of pyranometers not used for GTI estimation.
    """
    with open(file_path, 'r') as file:
        for _ in range(100):
            line = file.readline()
            if _ == line_number:
                match = re.search(r'\[(.*?)\]', line)
                if match:
                    content = match.group(1)
                    return [item.strip("'") for item in content.split(', ')]
    return None

def most_frequent_pyrano_combination(series):
    """
    Find the most frequent pyranometer combination in a series.

    Parameters:
        series (pd.Series): Series of pyranometer combinations.

    Returns:
        tuple: Most frequent combination and its percentage.
    """
    # Drop any NaN values in the Series
    non_empty_series = series.dropna()
    # Count the occurrences of each tuple
    value_counts = non_empty_series.value_counts()
    if value_counts.empty:
        return None, 0.0  # Return None or some default if the series is empty after dropping NaNs
    # Get the most frequent tuple and its count
    most_frequent_tuple = value_counts.idxmax()
    count = value_counts.max()
    # Calculate the percentage
    percentage = (count / len(non_empty_series)) * 100
    return most_frequent_tuple, round(percentage)

def create_gti_estimation_label(validation_pyrano):
    """
    Create a label for GTI estimation based on pyranometer names.

    Parameters:
        validation_pyrano (list): List of pyranometer names.

    Returns:
        dict: GTI estimation labels.
    """
    gti_estimation_label = {}
    direction_to_azimuth = {
        'N': 0,
        'NE': 45,
        'E': 90,
        'SE': 135,
        'S': 180,
        'SW': 225,
        'W': 270,
        'NW': 315
    }
    for item in validation_pyrano:
        direction, tilt = item.split('_')
        azimuth = direction_to_azimuth[direction]
        gti_estimation_label[item] = f'gti{azimuth}_{tilt}'
    return gti_estimation_label

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
