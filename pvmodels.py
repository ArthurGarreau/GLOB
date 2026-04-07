# -*- coding: utf-8 -*-
"""
Decomposition Models for Solar Radiation: Diffuse Fraction Calculation
======================================================================

This script provides implementations of several empirical models for estimating the diffuse fraction
of solar radiation from global horizontal irradiance (GHI) and clearness index (kt). These models are
widely used in solar energy applications to decompose global irradiance into direct and diffuse components.

Functions:
    - calculate_I0: Calculates extraterrestrial irradiance on a horizontal surface.
    - spencer: Spencer model for diffuse fraction estimation.
    - skartveit1: Skartveit-Olseth model for diffuse fraction estimation.
    - muneer2: Muneer model for diffuse fraction estimation.
    - mondol2: Mondol model for diffuse fraction estimation.
    - erbs: Erbs model for diffuse fraction estimation.
    - orgill_holland: Orgill and Hollands model for diffuse fraction estimation.


Created on Mon Feb  9 15:29:48 2026
@author: arthurg
"""

import numpy as np
import pandas as pd
import pvlib
from pvlib.location import Location
from pvlib.solarposition import get_solarposition

def calculate_I0(times, latitude, longitude):
    """
    Calculate extraterrestrial irradiance on a horizontal surface.

    Parameters:
        times (pd.DatetimeIndex): Time index for calculation.
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.

    Returns:
        pd.Series: Extraterrestrial irradiance (I0).
    """
    site = Location(latitude, longitude)
    solar_position = site.get_solarposition(times=times)
    dni_extra = pvlib.irradiance.get_extra_radiation(times)
    I_0 = dni_extra * np.cos(np.radians(solar_position['apparent_zenith']))
    return I_0

def spencer(kt, latitude_deg):
    """
    Spencer model for diffuse fraction estimation.

    Reference:
        Spencer, J. W. (1982). Solar Energy, 29(1), 19–32. https://doi.org/10.1016/0038-092X(82)90277-8
    """
    b = 0.940 + 0.0118 * np.abs(latitude_deg)
    c = 1.185 + 0.0135 * np.abs(latitude_deg)
    d = b - 0.75 * c
    a = b - 0.3 * c

    kd = np.select([kt <= 0.35, (kt > 0.35) & (kt <= 0.75), kt > 0.75], [a, b - c * kt, d])

    kd[kd > 1] = 1; kd[kd < 0] = 0

    if not isinstance(kd, pd.Series):
        kd = pd.Series(kd)
        kd.index = kt.index

    return kd

def skartveit1(kt, solarelevation_deg):
    """
    Skartveit-Olseth model for diffuse fraction estimation.

    Reference:
        Skartveit, Arvid, and Jan Asle Olseth. ‘Modelling Slope Irradiance at High Latitudes’. Solar Energy 36, no. 4 (1986): 333–44. https://doi.org/10.1016/0038-092X(86)90151-9


    """
    k1 = 0.87 - 0.56 * np.exp(-0.06 * solarelevation_deg)
    k2 = 1.09 * k1
    d1 = 0.15 + 0.43 * np.exp(-0.06 * solarelevation_deg)
    kc = 0.5 * (1 + np.sin((np.pi) * ((kt - 0.2) / (k1 - 0.2) - 0.5)))
    kc1 = 0.5 * (1 + np.sin((np.pi) * ((k2 - 0.2) / (k1 - 0.2) - 0.5)))
    a = 0.27
    b = 0

    cases = [
        (kt < 0.20),
        (0.20 <= kt) & (kt <= k2),
        (kt > k2)
    ]

    choices = [
        1,
        1 - (1 - d1) * (a * np.sqrt(kc) + b * kc + (1 - a - b) * kc ** 2),
        1 - k2 * (1 - (1 - (1 - d1) * (a * np.sqrt(kc1) + b * kc1 + (1 - a - b) * kc1 ** 2))) / kt
    ]

    kd = np.select(cases, choices)

    kd[kd > 1] = 1; kd[kd < 0] = 0

    if not isinstance(kd, pd.Series):
        kd = pd.Series(kd)
        kd.index = kt.index

    return kd

def muneer2(kt):
    """
    Muneer model for diffuse fraction estimation.

    Reference:
        Muneer, T., & Saluja, G. S. (1986). Building Services Engineering Research & Technology, 7(1), 37–43. https://doi.org/10.1177/014362448600700106
    """
    kd = np.select([(kt < 0.2), (0.2 <= kt)],
                   [0.98, 0.687 + 2.932 * kt - 8.546 * kt ** 2 + 5.227 * kt ** 3])

    kd[kd > 1] = 1; kd[kd < 0] = 0

    if not isinstance(kd, pd.Series):
        kd = pd.Series(kd)
        kd.index = kt.index

    return kd

def mondol2(kt):
    """
    Mondol model for diffuse fraction estimation.

    Reference:
        Mondol, J. D., Yohanis, Y. G., & Norton, B. (2008). Renewable Energy, 33(5), 1109–1120. https://doi.org/10.1016/j.renene.2007.06.005
    """
    kd = np.select(
        [(kt <= 0.20),
         (0.20 < kt) & (kt <= 0.7),
         (kt > 0.7)],
        [0.98,
         0.61092 + 3.6259 * kt - 10.171 * kt**2 + 6.338 * kt**3,
         0.672 - 0.474 * kt]
    )

    kd[kd > 1] = 1; kd[kd < 0] = 0

    if not isinstance(kd, pd.Series):
        kd = pd.Series(kd)
        kd.index = kt.index

    return kd

def erbs(kt):
    """
    Erbs model for diffuse fraction estimation.

    Reference:
        Erbs, D. G., Klein, S. A., & Duffie, J. A. (1982). Solar Energy, 28, 293–302. https://doi.org/10.1016/0038-092X(82)90302-4
    """
    kd = np.select([(kt <= 0.22), (0.22 < kt) & (kt <= 0.8), (kt > 0.8)],
                   [1.0 - 0.09 * kt, 0.9511 - 0.1604 * kt + 4.388 * kt**2 - 16.638 * kt**3 + 12.336 * kt**4,
                    0.165])

    kd[kd > 1] = 1; kd[kd < 0] = 0

    if not isinstance(kd, pd.Series):
        kd = pd.Series(kd)
        kd.index = kt.index

    return kd

def orgill_holland(kt):
    """
    Orgill and Hollands model for diffuse fraction estimation.

    Reference:
        Orgill, J. F., & Hollands, K. G. T. (1977). Solar Energy, 19, 357–359. https://doi.org/10.1016/0038-092X(77)90006-8
    """
    kd = np.select([(kt <= 0.35), (0.35 < kt) & (kt <= 0.75), (kt > 0.75)],
                   [1.0 - 0.249 * kt, 1.557 - 1.84 * kt, 0.177])

    kd[kd > 1] = 1; kd[kd < 0] = 0

    if not isinstance(kd, pd.Series):
        kd = pd.Series(kd)
        kd.index = kt.index

    return kd
