# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 17:26:46 2026

@author: arthurg
"""

import numpy as np

def get_mae_coefficient(azimuth_orientations, mae_data, azimuth, inclination, npyr):
    """Get mae coefficient for a given azimuth and inclination."""
    closest_inclination = min([45, 90, 135], key=lambda x: abs(x - inclination))
    if azimuth in ['NE', 'SE', 'SW', 'NW']:
        col = f"{azimuth}_{closest_inclination}"
        return mae_data.loc['nMAE_GLOB' + npyr, col] / 100
    elif azimuth in ['N', 'E', 'S', 'W']:
        aux = {'N': 'NE', 'E': 'SE', 'S': 'SW', 'W': 'NW'}
        col = f"{aux[azimuth]}_{closest_inclination}"
        return mae_data.loc['nMAE_GLOB25', col] / 100
    else:
        return np.nan

def calculate_irradiance_avg(df_estim, azimuth_angles, inclination_angles, idx):
    """Calculate irradiance_avg for monofacial or bifacial configurations."""
    irradiance_avg = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)
    if idx == 0:  # Monofacial
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                var_name = f"gti{azimuth}_{inclination}"
                if var_name in df_estim.columns:
                    irradiance_avg[i, j] = np.nanmean(df_estim[var_name])
    elif idx == 1:  # Bifacial
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                azimuth_ref = (azimuth + 180) % 360
                var_name = f"gti{azimuth}_{inclination}"
                var_name_ref = f"gti{azimuth_ref}_{180-inclination}"
                if var_name in df_estim.columns:
                    irradiance_avg[i, j] = np.nanmean(df_estim[var_name]) + np.nanmean(df_estim[var_name_ref])
    irradiance_avg[:, -1] = irradiance_avg[:, 0]  # Copy first to last
    return irradiance_avg

def apply_mae_adjustment(df_estim, azimuth_orientations, azimuth_angles, inclination_angles, mae_data, idx, npyr):
    """Apply nMAE adjustment to irradiance_avg."""
    irradiance_avg = calculate_irradiance_avg(df_estim, azimuth_angles, inclination_angles, idx)
    irradiance_avg_mae = np.zeros_like(irradiance_avg)
    for i, inclination in enumerate(inclination_angles):
        for j, azimuth in enumerate(azimuth_angles[:-1]):
            orientation = azimuth_orientations[j]
            if idx == 0:  # Monofacial
                mae_factor = get_mae_coefficient(azimuth_orientations, mae_data, orientation, inclination, npyr)
                irradiance_avg_mae[i, j] = irradiance_avg[i, j] * mae_factor
            elif idx == 1:  # Bifacial
                azimuth_ref = (azimuth + 180) % 360
                index = np.where(azimuth_ref == azimuth_angles)[0][0]
                orientation_ref = azimuth_orientations[index]
                mae_factor = get_mae_coefficient(azimuth_orientations, mae_data, orientation, inclination, npyr)
                mae_factor_ref = get_mae_coefficient(azimuth_orientations, mae_data, orientation_ref, 180 - inclination, npyr)
                irradiance_avg_mae[i, j] = (
                    np.nanmean(df_estim[f"gti{azimuth}_{inclination}"]) * mae_factor +
                    np.nanmean(df_estim[f"gti{azimuth_ref}_{180-inclination}"]) * mae_factor_ref
                )
        irradiance_avg_mae[i, -1] = irradiance_avg_mae[i, 0]
    return irradiance_avg_mae
