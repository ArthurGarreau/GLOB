
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
from scipy.optimize import least_squares
from joblib import Parallel, delayed
import re



def calculate_solar_angles(timestamps, latitude, longitude, altitude=6, temperature=-6):
    """
    Calculate solar angles (zenith, azimuth, declination, etc.) using pvlib.

    Parameters:
        timestamps (pd.DatetimeIndex): Timestamps for calculation.
        longitude (float): Longitude of the location.
        altitude (float): Altitude of the location (default: 6 m).
        temperature (float): Temperature (default: -6°C).

    Returns:
        pd.DataFrame: Solar position data (zenith, azimuth, declination, etc.).
    """
    
    solar_position = pvlib.solarposition.get_solarposition(
        time=timestamps, latitude=latitude, longitude=longitude, altitude=altitude, temperature=temperature)
                
    if isinstance(timestamps, pd.core.indexes.datetimes.DatetimeIndex):
        if len(timestamps)<=1:
            raise ValueError(f"Not the correct format of timestamps: {type(timestamps)}")
        else:
            eot = solar_position['equation_of_time']
            solar_position['hour_angle'] = pvlib.solarposition.hour_angle(
                timestamps, longitude, eot) # timestamps is already a pd.DatetimeIndex
    else:
        eot = solar_position['equation_of_time']
        solar_position['hour_angle'] = pvlib.solarposition.hour_angle(
            pd.DatetimeIndex([timestamps]), longitude, eot) # timestamps is not a pd.DatetimeIndex
                    
    solar_position['declination'] = np.degrees(pvlib.solarposition.declination_spencer71(timestamps.day_of_year))

    return solar_position

def incident_angle(solar_angles, plane_inclination, plane_azimuth, latitude, longitude):
    """
    Calculate the angle of incidence of the sun for a given plane.

    Parameters:
        solar_angles (pd.DataFrame): Dataframe of solar angles obtained with 
        the above function "calculate_solar_angles" (degrees).
        plane_inclination (float): Inclination angle of the plane (degrees).
        plane_azimuth (float): Azimuth angle of the plane (degrees).
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.

    Returns:
        tuple: Angle of incidence and zenith angle (degrees).
    """
    
    delta = np.radians(solar_angles['declination'].values); 
    omega = np.radians(solar_angles['hour_angle'].values); 
    beta = np.radians(plane_inclination)
    gamma = np.radians(plane_azimuth)
    phi = np.radians(latitude)
    
    theta_i = np.arccos(
        np.sin(delta) * np.sin(phi) * np.cos(beta) 
        - np.sin(delta) * np.cos(phi) * np.sin(beta) * np.cos(gamma) 
        + np.cos(delta) * np.cos(phi) * np.cos(beta) * np.cos(omega) 
        + np.cos(delta) * np.sin(phi) * np.sin(beta) * np.cos(omega) * np.cos(gamma) 
        + np.cos(delta) * np.sin(beta) * np.sin(omega) * np.sin(gamma)
   ) # see formula in Duffie, J. A., & Beckman, W. A. (2013). 
     # Available Solar Radiation Ch. 1. In Solar Engineering of Thermal Processes (pp. 43–137). 
     # John Wiley & Sons, Inc. https://doi.org/10.1002/9781118671603.ch2

    theta_i = np.degrees(theta_i)
    
    return theta_i
    
    
def Coef_b(solar_angles, plane_inclination, plane_azimuth, latitude, longitude):
    """
    Calculate the transposition coefficient for beam irradiance (b in function_least_square).

    Parameters:
        solar_angles (pd.DataFrame): Dataframe of solar angles obtained with 
        the function "calculate_solar_angles" (degrees).
        plane_inclination (float): Inclination angle of the plane (degrees).
        plane_azimuth (float): Azimuth angle of the plane (degrees).
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.

    Returns:
        float or array: Coef_b
    """

    theta_i = incident_angle(solar_angles, plane_inclination, plane_azimuth, latitude, longitude)
    theta_i[theta_i > 90] = 90; ####  !!!! VERY IMPORTANT !!!! The incident angles above 90 have to be set to 0 to avoid the beam to have an influence
    theta_i = np.radians(theta_i)
    theta_z = np.radians(solar_angles['zenith'].values)

    Coef_b = np.cos(theta_i) / np.cos(theta_z) 
    return Coef_b


def Coef_d(plane_inclination):
    """
    Calculate the transposition coefficient for diffuse irradiance (d in function_least_square).

    Parameters:
        plane_inclination (float): Inclination angle of the plane (degrees).
  
    Returns:
        float or array: Coef_d
    """
    beta = np.radians(plane_inclination)

    Coef_d = 0.5 * (1 + np.cos(beta))
    return Coef_d
    

def Coef_r(plane_inclination, plane_azimuth=None, solar_angles=None):
    """
    Calculate the transposition coefficient for reflected irradiance 
    (r_b and r_d in function_least_square).
    """
    beta = np.radians(plane_inclination)
    Coef_r = 0.5 * (1 - np.cos(beta))
    
    return Coef_r

def function_least_square(params, b, d, r, GTI, albedo=None):
    """
    Calculate the residuals for the least squares estimation of the quantity:
        GTI - (b + r*Albedo) * Beam_prime + (d + r*Albedo) * Diffuse_prime

    Parameters:
        params (list or array-like): The parameters to estimate 
            [Beam_prime, Diffuse_prime, Albedo] if albedo==None (this is the nonlinear case).
            [Beam_prime, Diffuse_prime]         if albedo is input (this is the linear case).
        b, d, r (array-like): The transpostion coefficients for Beam, Diffuse and Reflected.
        GTI (array-like): The observed values of GTI.

    Returns:
        array-like: The residuals.
    """
    if albedo==None: # nonlinear method case 
        x, y, z = params # (x,y,z) = (diffuse, beam, albedo)
        GTI_estim = (b + r * z) * y + (d + r * z) * x
    
    else: # linear method case
        x, y = params # (x,y) = (diffuse, beam)
        GTI_estim = (b + r * albedo) * y + (d + r * albedo) * x
        
    return GTI_estim - GTI
    


def estimation_diffuse_beam_Faiman(variables, glob_value, solar_angles, lat, lon, method='linear'):
    """
    Estimate diffuse and beam irradiance based on either linear or nonlinear equations.
    Parameters:
        variables (list): List of GLOB plane names.
        glob_value (pd.DataFrame): Global irradiance values.
        solar_angles (pd.DataFrame): Dataframe of solar angles.
        lat (float): Latitude of the location.
        lon (float): Longitude of the location.
        method (str): Method to use for estimation ('linear' or 'nonlinear').
    Returns:
        tuple: Estimated diffuse and beam irradiance, and error.
    """
    timestamps = glob_value.name
    albedo = glob_value['albedo']
    table_azim_incli = create_variable_table(variables)

    GTI_glob = glob_value[table_azim_incli.index].values

    # Check if GTI_glob is full of NaNs
    if np.all(np.isnan(GTI_glob)):
        return [np.nan, np.nan, np.nan, np.nan, np.nan]

    # Interpolate NaN values in GTI_glob
    GTI_glob = pd.Series(GTI_glob).interpolate(method='linear').values
    inclinations = table_azim_incli['inclination'].values
    azimuths = table_azim_incli['azimuth'].values

    b = np.array(Coef_b(solar_angles, inclinations, azimuths, lat, lon))
    d = np.array(Coef_d(inclinations))
    r = np.array(Coef_r(inclinations))

    X = np.column_stack((b, d, r))
    # Remove the NaN of the coefficients table
    rows_with_nan = np.isnan(X).any(axis=1)
    X = X[~rows_with_nan]
    b, d, r = X[:, 0], X[:, 1], X[:, 2]

    try:
        if method == 'linear':
            initial_guess = [10, 10]
            # Add bounds if applicable
            bounds = ([0, 0], [np.inf, np.inf])
            results = least_squares(
                function_least_square,
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
                function_least_square,
                initial_guess,
                bounds=bounds,
                args=(b, d, r, GTI_glob),
                max_nfev=1000  # Increase maximum number of function evaluations
            )
            D_prime, B_prime, Albedo = results.x

        # Check if the solution is valid
        if not results.success:
            return [np.nan, np.nan, np.nan, np.nan, np.nan]

    except Exception as e:
        print(f"An error occurred during optimization: {e}")
        return [np.nan, np.nan, np.nan, np.nan, np.nan]

    zenith = solar_angles['zenith'].values # degrees
    cos_z = np.cos(np.radians(zenith))
    I_0 = pvlib.irradiance.get_extra_radiation(timestamps)

    if not (np.isnan(D_prime) and np.isnan(B_prime)):
        D, B = solve_for_D_and_B(I_0, cos_z, D_prime, B_prime)

    D = D[0] if np.shape(D) == (1,) and not np.isnan(D[0]) else D
    B = B[0] if np.shape(B) == (1,) and not np.isnan(B[0]) else B
    D_prime = D_prime[0] if np.shape(D_prime) == (1,) and not np.isnan(D_prime[0]) else D_prime
    B_prime = B_prime[0] if np.shape(B_prime) == (1,) and not np.isnan(B_prime[0]) else B_prime
    albedo = albedo[0] if np.shape(albedo) == (1,) and not np.isnan(albedo[0]) else albedo

    if D > 0 and B > 0 and zenith < 90:
        return [np.round(D), np.round(B), np.round(albedo, 2), np.round(D_prime), np.round(B_prime)]
    else:
        return [np.nan, np.nan, np.nan, np.nan, np.nan]



def find_best_estimation(combs, glob_value, solar_angles, lat, lon, true_estimation, method='linear'):
    """
    Find the best combination of variables for the least-square estimation of 
    beam and diffuse (and albedo) using a parralelisation method.

    Parameters:
        combs (list): List of variable combinations.
        glob_value (pd.DataFrame): Global irradiance values.
        zenith_angle (float): Solar zenith angle (degrees).
        lat (float): Latitude of the location.
        lon (float): Longitude of the location.
        true_estimation (list of lists): List of measured [diffuse, beam].

    Returns:
        tuple: Best estimated diffuse and beam irradiance, and the best combination.
    """
    # Evaluate all combinations in parallel
    results = np.array( 
        Parallel(n_jobs=-1)(delayed(estimation_diffuse_beam_Faiman)(
        list(comb), glob_value, solar_angles, lat, lon, method
    ) for comb in combs) )
    
    true_diffuse, true_beam = true_estimation[0], true_estimation[1]
    # Extract errors from valid results
    errors = [np.sqrt( (results[:,0]-true_diffuse)**2 + (results[:,1]-true_beam)**2)]
    
    # Check if all results are NaN
    if np.all(np.isnan(errors)):
        return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan] # We add a Nan for the pyrano combination
    
    else:
        # Find the index of the minimum error among valid results
        best_index = np.nanargmin(errors)
        
        # Get the best result
        D, B, albedo, D_prime, B_prime = results[best_index,:]
        comb_opt = combs[best_index]
    
        return [np.round(D), np.round(B), np.round(albedo, 2), np.round(D_prime), np.round(B_prime), comb_opt]


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


def solve_for_D_and_B(I_0, cos_z, D_prime, B_prime):
    """
    Solve for D and B using the quadratic equation.
    See equation 1.11 (Appendix A) in Faiman et al., (1987).
    
    Parameters:
        I_0 (float): Extraterrestrial irradiance.
        cos_z (float): Cosine of the zenith angle.
        D_prime (float): Intermediate value for diffuse irradiance.
        I_prime (float): Intermediate value for beam irradiance.

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
        # print("The discriminant is negative, no real roots exist.")
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
            # print("No positive root found for D.")
            D, B = np.nan, np.nan
         
    return D, B

def calc_Dprime_Iprime(I_0, zenith, D, I):
    """
    Calculate D' and I' from D and I using the provided formulas.
    See equation 1.9 (Appendix A) in Faiman et al., (1987).

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

def azimuth_to_orientation(azimuth):
    """
    Convert an azimuth angle to its corresponding cardinal orientation.

    This function takes an azimuth angle as input and returns the cardinal orientation
    as a string based on predefined ranges.

    Parameters:
    azimuth (float): The azimuth angle in degrees, which can be positive or negative.

    Returns:
    str: The cardinal orientation corresponding to the azimuth angle. Possible return values
         are 'S', 'SW', 'W', 'NW', 'N', 'NE', 'E', 'SE', or 'Unknown' if the azimuth does not
         match any predefined range.
    """
    # Determine the orientation based on the azimuth
    if azimuth == 0:
        return 'S'
    elif azimuth == -45:
        return 'SE'
    elif azimuth == -90:
        return 'E'
    elif azimuth == -135:
        return 'NE'
    elif azimuth == 180:
        return 'N'
    elif azimuth == 135:
        return 'NW'
    elif azimuth == 90:
        return 'W'
    elif azimuth == 45:
        return 'SW'
    else:
        return 'Unknown'

def orientation_to_azimuth(orientation):
    """
    Convert a cardinal orientation to its corresponding azimuth angle.

    This function takes a cardinal orientation as input and returns the azimuth angle
    in degrees based on predefined values.

    Parameters:
    orientation (str): The cardinal orientation. Possible values are 
    'S', 'SW', 'W', 'NW', 'N', 'NE', 'E', and 'SE'.

    Returns:
    float: The azimuth angle in degrees corresponding to the cardinal orientation.
           Returns None if the orientation does not match any predefined value.
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
    
    reference_list = [
    'GHI', 'N_45', 'N_90', 'N_135', 'NE_45', 'NE_90', 'NE_135', 'E_45', 'E_90',
    'E_135', 'SE_45', 'SE_90', 'SE_135', 'S_45', 'S_90', 'S_135', 'SW_45',
    'SW_90', 'SW_135', 'W_45', 'W_90', 'W_135', 'NW_45', 'NW_90', 'NW_135']
    
    # Convert both lists to sets
    input_set = set(input_list)
    # Find the difference between the reference set and the input set
    missing_elements = [item for item in reference_list if item not in input_set]
    
    # Convert the result back to a list
    return missing_elements


def find_content_between_braces(file_path, line_number):
    """
    Function to return the name of the pyranometers that have not been used for
    the the estimations of GTI.
    ----------
    Returns
    list
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

