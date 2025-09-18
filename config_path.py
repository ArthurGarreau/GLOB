# -*- coding: utf-8 -*-
"""
Modify the path of the location of your scripts before executing.

Created on Tue May  6 11:52:03 2025
@author: arthurg
"""

from pathlib import Path


#### MAIN PATHS ####
SCRIPT_PATH = Path(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\GLOB_scripts")
# Location to output plots
PLOT_PATH = SCRIPT_PATH.parent 
# Location to retrieve data for estimations
DATA_PATH = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB")

#### SECONDARY PATHS ####
GLOB_DATA_PATH = DATA_PATH / "GLOB_Measurements"
NYA_DATA_PATH = DATA_PATH / "BSRN_NYA_Measurements"
B_D_DATA_PATH = DATA_PATH / "Beam_Diffuse_Estimations"
GTI_DATA_PATH = DATA_PATH / "GTI_Estimations"

HIGH_RES_PLOT_PATH = PLOT_PATH / "Figure_high_res_pdf"
LOW_RES_PLOT_PATH = PLOT_PATH / "Figure_low_res_png"
ALL_PLOT_PATH = PLOT_PATH / "Figure_all"

# Ensure directories exist
for path in [HIGH_RES_PLOT_PATH, LOW_RES_PLOT_PATH, ALL_PLOT_PATH, B_D_DATA_PATH, GTI_DATA_PATH]:
    path.mkdir(parents=True, exist_ok=True)