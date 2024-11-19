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


def calculate_hour_angle(time, latitude, longitude, timezone='UTC'):
    """
    Calculate the hour angle of the Sun for a given time, latitude, and longitude.

    Parameters:
    - time: pd.Timestamp or str, Time of observation (in UTC).
    - latitude: float, Latitude of the location.
    - longitude: float, Longitude of the location.
    - timezone: str, Time zone of the location (default is 'UTC').

    Returns:
    - hour_angle: float, Hour angle of the Sun in degrees.
    """
    # Convert the time to a pandas Timestamp if it's a string
    if isinstance(time, str):
        time = pd.to_datetime(time)

    # Create a timezone-aware datetime object
    time = time.tz_localize('UTC').tz_convert(timezone)
    
    # Convert date_str to a pandas Timestamp (at midnight)
    start_time = pd.to_datetime(time[0].date())
    
    # Generate a time range for every minute of the given day
    time_vector = pd.date_range(start=start_time, periods=1440, freq='min')  # 1440 minutes in a day
    solar_position = pvlib.solarposition.get_solarposition(time_vector, latitude, longitude)
    solar_noon = solar_position['zenith'].idxmin()
    solar_noon = solar_noon.tz_localize('UTC')


    # Calculate the hour angle (H = 15° * timedelta)
    hour_angle = (time - solar_noon) / pd.Timedelta(hours=1) * 15
    
    return hour_angle

def incidence_angle(phi, beta, delta, omega, gamma):
    """
    Calculate the cosine of the angle of incidence, θ.
    """
    theta = np.arccos(
        np.sin(delta) * np.sin(phi) * np.cos(beta) 
        - np.sin(delta) * np.cos(phi) * np.sin(beta) * np.cos(gamma) 
        + np.cos(delta) * np.cos(phi) * np.cos(beta) * np.cos(omega) 
        + np.cos(delta) * np.sin(phi) * np.sin(beta) * np.cos(omega) * np.cos(gamma) 
        + np.cos(delta) * np.sin(beta) * np.sin(omega) * np.sin(gamma)
    )
    return theta

def C_B_lambda(phi, beta, delta, omega, gamma):
    """
    Calculate the geometry coefficient for Beam irradiance (C_{B, lambda}).
    """
    cos_theta = np.array(np.cos(incidence_angle(phi, beta, delta, omega, gamma)))
    # cos_theta[cos_theta<0] = 0
    return cos_theta  # Only positive values contribute to irradiance

def C_D_lambda(phi, beta, delta, omega, gamma):
    """
    Calculate the geometry coefficient for Diffuse irradiance (C_{D, lambda}).
    Assumes isotropic diffuse model.
    """
    return 0.5 * (1 + np.cos(beta)) * (1 + np.sin(beta/2)) * (1 + C_B_lambda(phi, beta, delta, omega, gamma)**2 
                                                              * np.sin(incidence_angle(phi, 0, delta, omega, gamma)))

def C_R_lambda(beta):
    """
    Calculate the geometry coefficient for Reflected irradiance (C_{R, lambda}).
    Assumes isotropic ground reflection.
    """
    return 0.5 * (1 - np.cos(beta))


def I_B_n(df_T, C_B_df, C_D_df, C_R_df, rho):
    # Extract irradiance values
    I_T_E = df_T["E_45"]
    I_T_S = df_T["S_45"]
    I_T_W = df_T["W_45"]
    I_T_N = df_T["N_45"]
    I_T_h = df_T["GHI"]

    # Extract coefficients for each direction
    C_B_E = C_B_df["East"]
    C_B_S = C_B_df["South"]
    C_B_W = C_B_df["West"]
    C_B_N = C_B_df["North"]

    C_D_E = C_D_df["East"]
    C_D_S = C_D_df["South"]
    C_D_W = C_D_df["West"]
    C_D_N = C_D_df["North"]

    C_R_E = C_R_df["East"]
    C_R_S = C_R_df["South"]
    C_R_W = C_R_df["West"]
    C_R_N = C_R_df["North"]

    # Calculate numerators and denominators including North component
    numerator1 = ((I_T_E - I_T_h * rho * C_R_E) * C_D_S - (I_T_S - I_T_h * rho * C_R_S) * C_D_E)
    denominator1 = C_D_S * C_B_E - C_B_S * C_D_E
    
    numerator2 = ((I_T_S - I_T_h * rho * C_R_S) * C_D_W - (I_T_W - I_T_h * rho * C_R_W) * C_D_S)
    denominator2 = C_D_W * C_B_S - C_B_W * C_D_S
    
    numerator3 = ((I_T_W - I_T_h * rho * C_R_W) * C_D_N - (I_T_N - I_T_h * rho * C_R_N) * C_D_W)
    denominator3 = C_D_N * C_B_W - C_B_N * C_D_W
    
    numerator4 = ((I_T_N - I_T_h * rho * C_R_N) * C_D_E - (I_T_E - I_T_h * rho * C_R_E) * C_D_N)
    denominator4 = C_D_E * C_B_N - C_B_E * C_D_N
    

    # Initialize a result Series with the same index
    result = pd.Series(index=df_T.index, dtype=float)
    

    # Calculate values for different time ranges
    result.loc[(result.index.hour >= 11) & (result.index.hour < 17)] = (numerator1 + numerator2) / (denominator1 + denominator2)
    result.loc[(result.index.hour >= 17) & (result.index.hour < 23)] = (numerator2 + numerator3) / (denominator2 + denominator3)
    result.loc[(result.index.hour >= 5) & (result.index.hour < 11)] = (numerator3 + numerator4) / (denominator3 + denominator4)
    result.loc[(result.index.hour >= 23) | ((result.index.hour >= 0) & (result.index.hour < 5))] = (numerator4 + numerator1) / (denominator4 + denominator1)

    return result



def I_D(df_T, phi, beta, delta, omega, gamma):
    # Extract irradiance values
    I_T_E = df_T["E_90"]
    I_T_S = df_T["S_90"]
    I_T_W = df_T["W_90"]
    I_T_N = df_T["N_90"]

    C_D_df = pd.DataFrame({
        "East": C_D_lambda(phi, beta, delta, omega, gamma["East"]),
        "South": C_D_lambda(phi, beta, delta, omega, gamma["South"]),
        "West": C_D_lambda(phi, beta, delta, omega, gamma["West"]),
        "North": C_D_lambda(phi, beta, delta, omega, gamma["North"])
    }, index=df_T.index)

    # Initialize a result Series with the same index
    result = pd.Series(index=df_T.index, dtype=float)
    

    # Calculate values for different time ranges
    result.loc[(result.index.hour >= 20) | ((result.index.hour >= 0) & (result.index.hour < 2))] = I_T_S / C_D_df['South']
    result.loc[(result.index.hour >= 2) & (result.index.hour < 8)] = I_T_W / C_D_df['West']
    result.loc[(result.index.hour >= 8) & (result.index.hour < 14)] = I_T_N / C_D_df['North']
    result.loc[(result.index.hour >= 14) & (result.index.hour < 20)] = I_T_E / C_D_df['East']

    return result


def I_R(df_T, beta, rho):
    # Extract irradiance values
    I_T_h = df_T["GHI"]


    C_R_df = pd.DataFrame({
        "East": C_R_lambda(beta),
        "South": C_R_lambda(beta),
        "West": C_R_lambda(beta),
        "North": C_R_lambda(beta)},
    index=df_T.index)

    # Initialize a result Series with the same index
    result = pd.Series(index=df_T.index, dtype=float)
    

    # Calculate values for different time ranges
    result.loc[(result.index.hour >= 20) | ((result.index.hour >= 0) & (result.index.hour < 2))] = I_T_h * C_R_df['South'] * rho
    result.loc[(result.index.hour >= 2) & (result.index.hour < 8)] = I_T_h * C_R_df['West'] * rho
    result.loc[(result.index.hour >= 8) & (result.index.hour < 14)] = I_T_h * C_R_df['North'] * rho
    result.loc[(result.index.hour >= 14) & (result.index.hour < 20)] = I_T_h * C_R_df['East'] * rho

    return result



# def I_B_n_ES(I_T_E, I_T_S, I_T_h, C_B_E, C_D_E, C_R_E, C_B_S, C_D_S, C_R_S, rho):
#     """
#     Calculate I_{B,n}(ES), the beam normal irradiance between East and South.
#     """
#     numerator = ((I_T_E - I_T_h * rho * C_R_E) * C_D_S - (I_T_S - I_T_h * rho * C_R_S) * C_D_E)
#     denominator = C_D_S * C_B_E - C_B_S * C_D_E
    
#     return numerator / denominator 

# def I_B_n_SW(I_T_S, I_T_W, I_T_h, C_B_S, C_D_S, C_R_S, C_B_W, C_D_W, C_R_W, rho):
#     """
#     Calculate I_{B,n}(SW), the beam normal irradiance between South and West.
#     """
#     numerator = ((I_T_S - I_T_h * rho * C_R_S) * C_D_W - (I_T_W - I_T_h * rho * C_R_W) * C_D_S)
#     denominator = C_D_W * C_B_S - C_B_W * C_D_S
#     return numerator / denominator

# def I_B_n_WN(I_T_W, I_T_N, I_T_h, C_B_W, C_D_W, C_R_W, C_B_N, C_D_N, C_R_N, rho):
#     """
#     Calculate I_{B,n}(WN), the beam normal irradiance between West and North.
#     """
#     numerator = ((I_T_W - I_T_h * rho * C_R_W) * C_D_N - (I_T_N - I_T_h * rho * C_R_N) * C_D_W)
#     denominator = C_D_N * C_B_W - C_B_N * C_D_W
#     return numerator / denominator 

# def I_B_n_NE(I_T_N, I_T_E, I_T_h, C_B_N, C_D_N, C_R_N, C_B_E, C_D_E, C_R_E, rho):
#     """
#     Calculate I_{B,n}(NE), the beam normal irradiance between North and East.
#     """
#     numerator = ((I_T_N - I_T_h * rho * C_R_N) * C_D_E - (I_T_E - I_T_h * rho * C_R_E) * C_D_N)
#     denominator = C_D_E * C_B_N - C_B_E * C_D_N
#     return numerator / denominator 

