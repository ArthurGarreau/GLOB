# -*- coding: utf-8 -*-
"""

Created on Tue May  6 11:52:03 2025
@author: arthurg
"""

from pathlib import Path

#### MAIN PATHS ####
SCRIPT_PATH = Path(".")
# Location to output plots
FIG_PATH = SCRIPT_PATH / ".." / "FIGURES" 
# Location to retrieve data for estimations
DATA_PATH = SCRIPT_PATH / "DATA_toshare" 

#### SECONDARY PATHS ####
B_D_DATA_PATH = DATA_PATH / "Estim_Beam_Diffuse"
GTI_DATA_PATH = DATA_PATH / "Estim_GTI"

HIGH_RES_PLOT_PATH = FIG_PATH / "Figure_high_res"
LOW_RES_PLOT_PATH = FIG_PATH / "Figure_low_res"
ALL_PLOT_PATH = FIG_PATH / "Figure_all"

# Ensure directories exist
for path in [HIGH_RES_PLOT_PATH, LOW_RES_PLOT_PATH, ALL_PLOT_PATH,\
             B_D_DATA_PATH, GTI_DATA_PATH]:
    path.mkdir(parents=True, exist_ok=True)
    

################################ END ##########################################