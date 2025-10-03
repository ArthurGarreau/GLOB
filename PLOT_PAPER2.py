# -*- coding: utf-8 -*-
"""
Solar Irradiance Analysis and Plot Script
==================================================
This script performs a comprehensive analysis and visualization of solar irradiance components
using data from GLOB pyranometers and BSRN reference measurements at Ny-Ålesund.

Key Features:
-------------
- Compares beam and diffuse irradiance from BSRN and GLOB instruments
- Generates polar heatmaps of GTI (Global Tilted Irradiance) for different orientations and months
- Visualizes annual averages of monofacial and bifacial GTI
- Analyzes the most used pyranometer combinations for 3-sensor configurations
- Produces high-quality figures for publication

Figures Generated:
-----------------
1. Figure 4: Daily evolution of beam and diffuse components
2. Figures 5-6: Monthly average polar heatmaps of incoming and reflected GTI
3. Figure 7: Annual average GTI for monofacial and bifacial configurations (estimated)
4. Figure 8: Annual average GTI for monofacial and bifacial configurations (measured)
5. Figure 10: Heatmaps of most used pyranometer combinations by time of day

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: April 25, 2025
"""

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import calendar
import xarray as xr
import pvlib
import glob_functions_calculation as fct
from config_path import GTI_DATA_PATH, B_D_DATA_PATH, HIGH_RES_PLOT_PATH, \
LOW_RES_PLOT_PATH, ALL_PLOT_PATH, GLOB_DATA_PATH, NYA_DATA_PATH

# --- Global Parameters ---
method = 'nonlinear'  # Default method for irradiance decomposition
f = 10  # Data frequency in minutes
pyrano_nr = 5  # Number of pyranometers used in estimations

###############################################################################
# --- File Paths ---
# Define paths to input data files
gti_estimation_datafile = GTI_DATA_PATH / \
    f"2023-24_estimation_GTI_{f}min_{method}_{pyrano_nr}pyrano.csv"
gti_nya_datafile = GTI_DATA_PATH / \
    f"2025_estimation_GTI_{f}min_NYA.csv"
B_D_estimations_datafile = B_D_DATA_PATH / \
    f"2025_estimation_beam_diffuse_{f}min_{method}_{pyrano_nr}pyrano.csv"
bsrn_datafile = NYA_DATA_PATH / "NYA_radiation_2025-all.tab"
glob_data_file = GLOB_DATA_PATH / "GLOB_data_10min_2023-24.nc"

# %% Fig 5: Beam and diffuse daily evolution
# -----------------------------------------
# This section plots the daily evolution of beam and diffuse irradiance components
# for a specific date, comparing BSRN reference data with GLOB estimations and model predictions

# --- Parameters ---
date_str = '2025-03-29'  # Date for analysis
pyrano_nr = 25  # Number of pyranometers
method = 'linear'  # Method used for this specific figure
B_D_estimations_datafile = B_D_DATA_PATH / \
    f"2025_estimation_beam_diffuse_{f}min_{method}_{pyrano_nr}pyrano.csv"

# --- Helper Functions ---
def round_up_to_hundred(x):
    """Round up to the nearest hundred."""
    return np.ceil(x / 100) * 100

def calculate_decomposition_models(ghi, time, latitude, longitude):
    """Calculate beam and diffuse components using decomposition models."""
    solar_position = pvlib.solarposition.get_solarposition(time, latitude, longitude)
    zenith = solar_position['zenith'].values
    # Erbs model
    erbs = pvlib.irradiance.erbs(ghi, zenith, time)
    beam_erbs, diffuse_erbs = erbs['dni'], erbs['dhi']
    # Perez model
    perez_dni = pvlib.irradiance.dirint(ghi, zenith, time)
    beam_perez = perez_dni
    diffuse_perez = ghi - np.cos(np.radians(zenith)) * perez_dni
    # Orgill and Holland model
    orgill_holland = pvlib.irradiance.orgill_hollands(ghi, zenith, time)
    beam_oh, diffuse_oh = orgill_holland['dni'], orgill_holland['dhi']
    return beam_erbs, diffuse_erbs, beam_perez, diffuse_perez, beam_oh, diffuse_oh

# --- Load Data ---
# Load BSRN data
bsrn_data = pd.read_csv(bsrn_datafile, sep='\t',skiprows=24, parse_dates=['Date/Time'],
    index_col='Date/Time').resample(f'{f}min').first()
# Load GLOB data
glob_estim_data = pd.read_csv(B_D_estimations_datafile,parse_dates=['Timestamp'],
    index_col='Timestamp', sep='\t', header=10)
# Load GHI data from NetCDF file
with xr.open_dataset(GLOB_DATA_PATH / f"GLOB_data_10min_2025.nc") as ds:
    ghi_day = ds['GHI'].to_dataframe()
    latitude, longitude = ds.latitude.values, ds.longitude.values

# --- Filter Data for Specific Date ---
bsrn_data = bsrn_data.loc[date_str]
glob_data = glob_estim_data.loc[date_str]
ghi_day_current = ghi_day.loc[date_str]

# --- Extract Components ---
beam_bsrn = bsrn_data['DIR']
diffuse_bsrn = bsrn_data['DIF']
beam_glob = glob_data['Beam']
diffuse_glob = glob_data['Diffuse']

# --- Calculate Decomposition Models ---
time = pd.date_range(start=f"{date_str} 00:00:00", end=f"{date_str} 23:50:00",
    freq=f"{f}min")
beam_erbs, diffuse_erbs, beam_perez, diffuse_perez, beam_oh, diffuse_oh = calculate_decomposition_models(
    ghi_day_current['GHI'].values, time, latitude, longitude)

# --- Plot Data ---
plt.rc('font', size=14)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
# Beam subplot (top)
ax1.plot(beam_bsrn.index, beam_bsrn, 'k', marker='.', linewidth=2, label='BSRN Ny-Ålesund (truth)')
ax1.plot(beam_glob.index, beam_glob, 'royalblue', linewidth=2, label=f'Estimation with GLOB\n{pyrano_nr} pyranometers {method}')
ax1.plot(time, beam_perez, color='lightcoral', linestyle='--', label='Model Perez')
ax1.plot(time, beam_erbs, color='firebrick', linestyle='--', label='Model Erbs')
ax1.plot(time, beam_oh, color='maroon', linestyle='--', label='Model Orgill and Hollands')
ax1.set_title('Beam', fontsize=16, fontweight='bold')
ax1.set_ylabel('Irradiance [$W \ m^{-2}$]')
ax1.legend(fontsize=12, edgecolor='none')
# ax1.set_ylim(0, round_up_to_hundred(np.nanmax([beam_bsrn, beam_glob, beam_erbs, beam_oh, beam_perez])))
ax1.grid(True, linestyle=':')

# Diffuse subplot (bottom)
ax2.plot(diffuse_bsrn.index, diffuse_bsrn, 'k', marker='.', label='BSRN Ny-Ålesund (truth)')
ax2.plot(diffuse_glob.index, diffuse_glob, 'royalblue', linewidth=2, label=f'Estimation with GLOB\n{pyrano_nr} pyranometers {method}')
ax2.plot(time, diffuse_perez, color='lightcoral', linestyle='--', label='Model Perez')
ax2.plot(time, diffuse_erbs, color='firebrick', linestyle='--', label='Model Erbs')
ax2.plot(time, diffuse_oh, color='maroon', linestyle='--', label='Model Orgill and Hollands')
ax2.set_title('Diffuse', fontsize=16, fontweight='bold')
ax2.set_xlabel('Time (UTC)')
ax2.set_ylabel('Irradiance [$W \ m^{-2}$]')
# ax2.set_ylim(0, round_up_to_hundred(np.nanmax([diffuse_bsrn, diffuse_glob, diffuse_erbs, diffuse_oh, diffuse_perez])))
ax2.grid(True, linestyle=':')

# Format x-axis
start_time = pd.Timestamp(f"{date_str} 03:00:00")
end_time = pd.Timestamp(f"{date_str} 21:00:00")
ax1.set_xlim(start_time, end_time)
ax2.set_xlim(start_time, end_time)
ax2.xaxis.set_major_locator(mdates.HourLocator(interval=1))
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
ax2.annotate(date_str, xy=(0.90, -0.15), xycoords='axes fraction')
plt.tight_layout()

# Save the plot
# plt.savefig(ALL_PLOT_PATH / f"Beam_Diffuse_{date_str}_{pyrano_nr}pyrano_{method}.png.png", dpi=300, bbox_inches='tight')
plt.savefig(LOW_RES_PLOT_PATH / "Figure 5.png", dpi=300, bbox_inches='tight')
# plt.savefig(HIGH_RES_PLOT_PATH / "Figure 5.pdf", format='pdf', bbox_inches='tight')
# Optionally show the plot
plt.close()

# %% Fig 6-7: Combined polar heatmap of GTI - Monthly average of incoming and reflected GTI
# --------------------------------------------------------------------------------
# This section creates polar heatmaps showing monthly averages of GTI for sky-facing and ground-facing planes
# for the period from March to September in 2023-2024
# Moreover, this code enables saving the data from the plot in a .csv file.

# --- Load GTI Data ---
gti_data = pd.read_csv(gti_estimation_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)
gti_data = gti_data[~gti_data.index.duplicated(keep='first')]

# Define the monthly date ranges for 2023 and 2024
year = "2023-24"
monthly_date_ranges = [
    pd.date_range(
        start=f'2023-{month:02d}-01',
        end=f'2023-{month:02d}-{pd.Timestamp(f"2023-{month:02d}-01").days_in_month}',
        freq='10min'
    ).union(
        pd.date_range(
            start=f'2024-{month:02d}-01',
            end=f'2024-{month:02d}-{pd.Timestamp(f"2024-{month:02d}-01").days_in_month}',
            freq='10min'
        )
    )
    for month in range(4, 10)  # April to September (4=April, 9=September)
]

# Define the azimuth orientations and their corresponding angles in degrees
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
azimuth_angles = np.array([0, 45, 90, 135, 180, 225, 270, 315, 360])

# Define the inclination angles and titles
titles = [f"Sky-facing planes ({year})", f"Ground-facing planes ({year})"]
names = ["skyfacing", "groundfacing"]

# Initialize a list to store monthly irradiance matrices for CSV export
monthly_matrices = {name: [] for name in names}

for idx, title in enumerate(titles):
    if idx == 0:
        inclination_angles = np.arange(0, 90, 5)  # incoming solar light
    if idx == 1:
        inclination_angles = np.arange(90, 185, 5)  # reflected solar light

    # Initialize a grid to accumulate irradiance values for monthly averages
    theta, r = np.meshgrid(np.deg2rad(azimuth_angles), inclination_angles)

    # Create a 3x3 grid for the plots
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 4, figure=fig, width_ratios=[1, 1, 1, 0.1])  # Allocate space for the colorbar
    subgrid_count = [0, 1, 2, 4, 5, 6]

    # Initialize a list to store the contourf objects
    contourf_objects = []

    # Loop through each month's date range
    for idx_month, monthly_date_range in enumerate(monthly_date_ranges):
        # Determine the subplot position
        ax = fig.add_subplot(gs[subgrid_count[idx_month]], projection='polar')
        monthly_date_range = monthly_date_range.tz_localize('UTC')
        df_estim = gti_data.reindex(monthly_date_range)

        # Initialize a grid to accumulate irradiance values for the monthly average
        irradiance_avg = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)

        # Get the current month
        current_month = monthly_date_range[0].month

        # Populate the grid
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                var_name = f"gti{azimuth}_{inclination}"
                if var_name in df_estim.columns:
                    irradiance_avg[i, j] = np.nanmean(df_estim[var_name])

        # Duplicate first column to last for continuity
        irradiance_avg[:, -1] = irradiance_avg[:, 0]

        # Store the irradiance matrix for CSV export
        irradiance_df = pd.DataFrame(irradiance_avg, index=inclination_angles, columns=azimuth_orientations)
        irradiance_df.index.name = 'Tilt'
        irradiance_df.columns.name = 'Azimuth'
        irradiance_df = irradiance_df.stack().reset_index()
        irradiance_df['Month'] = calendar.month_name[current_month]
        irradiance_df['Case'] = names[idx]
        irradiance_df = irradiance_df.rename(columns={0: 'Irradiance'})
        monthly_matrices[names[idx]].append(irradiance_df)

        # Plot with improved color scaling
        vmax = np.nanmax(irradiance_avg)
        contour = ax.contourf(theta, r, irradiance_avg, cmap='gnuplot2', levels=np.arange(0, 340, 10), alpha=1)
        contourf_objects.append(contour)
        contour_lines = ax.contour(theta, r, irradiance_avg, colors='white',
                                   linewidths=0.8, levels=[100, 150, 200, 240, 250, 260], zorder=1)
        ax.clabel(contour_lines, inline_spacing=0.1, fontsize=9, fmt='%1.0f')

        # Add colorbar and formatting
        ax.set_xticks(np.radians(azimuth_angles))
        ax.set_xticklabels(azimuth_orientations, fontsize=12)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_yticks(inclination_angles[::2])
        ax.set_yticklabels(inclination_angles[::2], fontsize=12)
        ax.set_rlabel_position(150)
        ax.spines['polar'].set_color('black')
        ax.tick_params(axis='both', colors='black')
        ax.text(np.radians(150), inclination_angles[-1]+10, "Tilt angle (°)")
        ax.set_title(f'{calendar.month_name[current_month]}', fontsize=12, fontweight="bold")

    # Add a single colorbar to the figure
    vmin = min([co.get_array().min() for co in contourf_objects])
    vmax = max([co.get_array().max() for co in contourf_objects])
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cbar_ax = fig.add_subplot(gs[0:2, 3])
    cbar = fig.colorbar(contour, cax=cbar_ax)
    cbar.set_ticks(np.arange(0, 360, 30))
    fig.text(0.02, 0.96, title, fontsize=16, fontweight="bold", bbox=dict(facecolor='k', alpha=0.2))
    cbar.set_label(label='Irradiance [$W \ m^{-2}$]', size=14)
    cbar.ax.tick_params(labelsize=12)
    plt.tight_layout(rect=[0, 0, 0.9, 0.95])

    # Save the combined plot (uncomment if needed)
    # fig.savefig(ALL_PLOT_PATH / f"monthly_avg_polar_heatmap_{names[idx]}_{year}.png", dpi=300, bbox_inches='tight')
    # fig.savefig(LOW_RES_PLOT_PATH / f"Figure {6+idx}.png", dpi=300, bbox_inches='tight')
    # fig.savefig(HIGH_RES_PLOT_PATH / f"Figure {6+idx}.pdf", format='pdf', bbox_inches='tight')
    # plt.close()

# Combine all monthly matrices into a single DataFrame and save to CSV
all_months_df = pd.concat([pd.concat(monthly_matrices['skyfacing']), pd.concat(monthly_matrices['groundfacing'])])
all_months_pivot = all_months_df.pivot_table(
    index=['Month', 'Tilt', 'Case'],
    columns='Azimuth',
    values='Irradiance'
).reset_index()

# Round the irradiance values to 2 decimal places
all_months_pivot = all_months_pivot.round()
# Save to CSV
all_months_pivot.to_csv(HIGH_RES_PLOT_PATH / 'Figure 6-7.csv', sep='\t',index=False)

# %% Fig 8-9: Average GTI in 2023-24 on monofacial and bifacial
# -----------------------------------------------------------
# This section creates polar heatmaps showing annual averages of GTI for monofacial and bifacial configurations

# --- Load Data ---
glob_data_file = GLOB_DATA_PATH / f"GLOB_data_{f}min_2023-24.nc"
gti_data = pd.read_csv(gti_estimation_datafile, sep='\t', parse_dates=True, index_col='Timestamp', header=10)

# Create the date ranges
dates_2023 = pd.date_range('2023-04-01', '2023-09-30', freq='10min')
dates_2024 = pd.date_range('2024-04-01', '2024-09-30', freq='10min')
date_range = dates_2023.union(dates_2024); date_range = date_range.tz_localize('UTC')
gti_data = gti_data[~gti_data.index.duplicated(keep='first')]
df_estim = gti_data.reindex(date_range)

# Define the azimuth orientations and their corresponding angles in degrees
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
azimuth_angles = np.array([0, 45, 90, 135, 180, 225, 270, 315, 360])
# Define the inclination angles
inclination_angles = np.arange(0, 95, 5)  # Different sets of inclination angles

titles = ['(a) Monofacial\n     sky-facing\n     planes', '(b) Bifacial']
fig = plt.figure(figsize=(12,6))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 0.05])  # Allocate space for the colorbar

# Loop through each set of inclination angles
for idx, title in enumerate(titles):
    # Initialize a grid to accumulate irradiance values for the annual average
    irradiance_avg = np.zeros((len(inclination_angles), len(azimuth_angles)), dtype=float)
    # Populate the grid
    if idx == 0:
        levels = [150,200,210,215]
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                var_name = f"gti{azimuth}_{inclination}"
                if var_name in df_estim.columns:
                    irradiance_avg[i, j] = np.nanmean(df_estim[var_name])
    elif idx == 1:
        levels = [250,300,315, 317, 320]
        for i, inclination in enumerate(inclination_angles):
            for j, azimuth in enumerate(azimuth_angles):
                azimuth_ref = (azimuth + 180) % 360
                var_name = f"gti{azimuth}_{inclination}"
                var_name_ref = f"gti{azimuth_ref}_{180-inclination}"
                if var_name in df_estim.columns:
                    irradiance_avg[i, j] = np.nanmean(df_estim[var_name]) + np.nanmean(df_estim[var_name_ref])

    irradiance_avg[:,-1] = irradiance_avg[:,0]
     # Create the annual polar heatmap for the current set of inclination angles
    theta, r = np.meshgrid(np.deg2rad(azimuth_angles), inclination_angles)

    ax = fig.add_subplot(gs[idx], projection='polar')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    # Plot with improved color scaling
    contour = ax.contourf(theta, r, irradiance_avg, cmap='gnuplot2', levels=np.arange(0, 395, 5), alpha=1)
    contour_lines = ax.contour(theta, r, irradiance_avg, colors='white', linewidths=0.8, levels=levels, zorder=1)
    ax.clabel(contour_lines, inline_spacing=0.1, fontsize=10, manual=True, fmt='%1.0f')
    ax.set_xticks(np.radians(azimuth_angles))
    ax.set_xticklabels(azimuth_orientations, fontsize=12)
    ax.set_yticks(inclination_angles[0::2])
    ax.set_yticklabels(inclination_angles[0::2], fontsize=12)
    ax.set_rlabel_position(150)
    ax.spines['polar'].set_color('black')
    ax.tick_params(axis='both', colors='black')
    ax.text(np.radians(150), inclination_angles[-1] + 10, "Tilt angle [°]")
    ax.set_title(title, fontsize = 12, fontweight= "bold", loc='left', bbox=dict(facecolor='k', alpha=0.2))

cbar_ax = fig.add_subplot(gs[2])
cbar = fig.colorbar(contour, cax=cbar_ax)
cbar.set_ticks(np.arange(0, 420, 30))
fig.text(0.05, 0.95, "GLOB 5 pyranometers nonlinear estimations (2023-2024)", fontsize = 16, fontweight= "bold")
cbar.set_label(label='Irradiance [$W \ m^{-2}$]', size=14)  # Adjust the size as needed
cbar.ax.tick_params(labelsize=12)  # Adjust the size as needed
plt.tight_layout(rect=[0, 0, 0.9, 0.95])

# Save and show the monthly plot
fig.savefig(ALL_PLOT_PATH / "annual_avg_polar_heatmap_GLOB_estim_2023-24.png", dpi=300, bbox_inches='tight')
fig.savefig(LOW_RES_PLOT_PATH / "Figure 8.png", dpi=300, bbox_inches='tight')
fig.savefig(HIGH_RES_PLOT_PATH / "Figure 8.pdf", format='pdf', bbox_inches='tight')
# plt.close()

# %% Fig 9: GTI measured with GLOB (must run Fig 8 before)
# -------------------------------
# This section creates polar heatmaps showing annual averages of GTI based on actual GLOB measurements

# Create the date ranges
dates_2023 = pd.date_range('2023-04-01', '2023-09-30', freq='10min')
dates_2024 = pd.date_range('2024-04-01', '2024-09-30', freq='10min')
date_range = dates_2023.union(dates_2024); date_range = date_range.tz_localize('UTC')

# --- Load Data ---
ds_glob = xr.open_dataset(glob_data_file)
df_glob = ds_glob.to_dataframe(); df_glob.index = df_glob.index.tz_localize('UTC')
df_glob = df_glob.reindex(date_range, fill_value=pd.NA)
# We mask the glob values with the same NaN mask as in the estimations data
# to perform the averaging on the same data.
# mask = df_estim['gti0_45'].isna()
# df_glob = df_glob[~mask]

# Define azimuth orientations and angles
azimuth_orientations = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW','N']
azimuth_angles = np.array([0, 45, 90, 135, 180, 225, 270, 315, 360])
# Define inclination angles
inclination_angles = np.array([0, 45, 90])
# Initialize a grid to accumulate irradiance values
theta, r = np.meshgrid(np.deg2rad(azimuth_angles), inclination_angles)
irradiance_avg = np.zeros_like(theta, dtype=float)
titles = ['(a) Monofacial\n     sky-facing\n     planes', '(b) Bifacial']
fig = plt.figure(figsize=(12,6))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 0.05])  # Allocate space for the colorbar

for idx, title in enumerate(titles):
    # Populate the grid with average irradiance values
    for i, inclination in enumerate(inclination_angles):
        for j, azimuth in enumerate(azimuth_angles[:-1]):
            orientation = azimuth_orientations[j]
            orientation_ref = fct.azimuth_to_orientation((azimuth + 180) % 360)

            if idx == 0:
                levels = [150,200,210,215]
                if inclination == 0:
                    var_name = 'GHI'
                    var_name_ref = 'GHI_ground'
                    irradiance_avg[i, j] = np.nanmean(df_glob[var_name])

                else:
                    var_name = f"{orientation}_{inclination}"
                    irradiance_avg[i, j] = np.nanmean(df_glob[var_name])
            if idx == 1:
                levels = [200,250,300,305]
                if inclination == 0:
                    var_name = 'GHI'
                    var_name_ref = 'GHI_ground'
                    irradiance_avg[i, j] = np.nanmean(df_glob[var_name]) + np.nanmean(df_glob['albedo']*df_glob[var_name])
                else:
                    var_name = f"{orientation}_{inclination}"
                    var_name_ref = f"{orientation_ref}_{180-inclination}" # opposite plane
                    irradiance_avg[i, j] = np.nanmean(df_glob[var_name]) + np.nanmean(df_glob[var_name_ref])

    irradiance_avg[:,-1] = irradiance_avg[:,0]
    # Create the polar heatmap
    ax = fig.add_subplot(gs[idx], projection='polar')

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    # Plot with improved color scaling
    contour = ax.contourf(theta, r, irradiance_avg, cmap='gnuplot2', levels=np.arange(0, 395, 5), alpha=1)
    contour_lines = ax.contour(theta, r, irradiance_avg, colors='white', linewidths=0.8, levels=levels, zorder=1)
    ax.clabel(contour_lines, inline_spacing=0, manual=True, fontsize=10, fmt='%1.0f')
    ax.set_xticks(np.radians(azimuth_angles))
    ax.set_xticklabels(azimuth_orientations, fontsize=12)
    ax.set_yticks(inclination_angles)
    ax.set_yticklabels(inclination_angles, fontsize=12)
    ax.set_rlabel_position(150)
    ax.spines['polar'].set_color('black')
    ax.tick_params(axis='both', colors='black')
    ax.text(np.radians(150), inclination_angles[-1]+10, "Tilt angle (°)")
    ax.set_title(title, fontsize = 12, fontweight= "bold", loc='left', bbox=dict(facecolor='k', alpha=0.2))

cbar_ax = fig.add_subplot(gs[2])
cbar = fig.colorbar(contour, cax=cbar_ax)
cbar.set_ticks(np.arange(0, 420, 30))
fig.text(0.05, 0.95, "GLOB measurements (2023-2024)", fontsize = 16, fontweight= "bold")
cbar.set_label(label='Irradiance [$W \ m^{-2}$]', size=14)  # Adjust the size as needed
cbar.ax.tick_params(labelsize=12)  # Adjust the size as needed
plt.tight_layout(rect=[0, 0, 0.9, 0.95])

# Save the plot
# fig.savefig(ALL_PLOT_PATH / "annual_avg_polar_heatmap_GLOB_meas_2023-24.png", dpi=300, bbox_inches='tight')
fig.savefig(LOW_RES_PLOT_PATH / "Figure 9.png", dpi=300, bbox_inches='tight')
# fig.savefig(HIGH_RES_PLOT_PATH / "Figure 9.pdf", format='pdf', bbox_inches='tight')
# plt.close()

# %% Fig 10: Most used pyranometers for estimations with a combination of 3
# ------------------------------------------------------------------------
# This section creates heatmaps showing the most used pyranometer combinations
# for 3-sensor configurations throughout the day

# --- Load Data ---
beam_diffuse_datafile = B_D_DATA_PATH / "2025_bestestimation_beam_diffuse_10min_linear.csv"
data = pd.read_csv(beam_diffuse_datafile, parse_dates=['Timestamp'], index_col='Timestamp', sep='\t', header=10)

# Define the cardinal directions and tilt angles
cardinal_directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'GHI']
tilt_angles = [45, 90]

# Create a figure with a custom layout
fig = plt.figure(figsize=(12, 14))
plt.rcParams['font.size'] = 12
gs = GridSpec(5, 3, figure=fig, width_ratios=[1, 1, 0.1])  # Allocate space for the colorbar
subgrid_count = [0,3,6,9,1,4,7,10]

# Loop through each 3-hour period
for i in range(8):
    start_hour = i * 3
    end_hour = (i * 3 + 3) % 24
    start_time = f"{start_hour:02d}:00"
    end_time = f"{end_hour:02d}:00"
    # Create subplot for each period
    ax = fig.add_subplot(gs[subgrid_count[i]])
    # Filter the data between the specified hours
    filtered_data = data.between_time(start_time, end_time)
    # Extract the Pyrano_Combination column
    combinations = filtered_data['Pyrano_Combination']
    # Initialize a dictionary to count the occurrences of each combination
    combination_counts = {}
    # Count the occurrences of each combination
    for combo_list in combinations:
        if isinstance(combo_list, str):
            combo_list = combo_list.replace("'", "").strip("()")
            combos = combo_list.split(", ")
            for combo in combos:
                if combo in combination_counts:
                    combination_counts[combo] += 1
                else:
                    combination_counts[combo] = 1
    # Create a matrix to represent the counts for the heatmap
    heatmap_data = np.zeros((len(tilt_angles), len(cardinal_directions)))
    # Populate the heatmap data matrix
    for combo, count in combination_counts.items():
        if combo != 'GHI':
            direction, angle = combo.split('_')
        else:
            direction, angle = 'GHI', 90
        angle = int(angle)
        if direction in cardinal_directions and angle in tilt_angles:
            dir_index = cardinal_directions.index(direction)
            angle_index = tilt_angles.index(angle)
            heatmap_data[angle_index, dir_index] = count
    # Normalize the heatmap data
    if np.sum(heatmap_data) > 0:
        heatmap_data = heatmap_data / np.sum(heatmap_data) * 100
    heatmap_data[0,-1]= np.nan
    # Plot the heatmap
    c = ax.imshow(heatmap_data, cmap='viridis', aspect='auto', vmin=0, vmax=10)
    # Set the labels for the cardinal directions and tilt angles
    ax.set_xticks(np.arange(len(cardinal_directions)))
    ax.set_xticklabels(cardinal_directions)
    ax.set_yticks(np.arange(len(tilt_angles)))
    ax.set_yticklabels(tilt_angles)
    # Set the title for the subplot
    ax.set_title(f'({chr(97 + i)}) {start_time} - {end_time} UTC', fontsize=14, fontweight='bold')
    if i < 4: ax.set_ylabel('Tilt angle [°]')

# Create a larger subplot for the full day at the bottom spanning both columns
ax_full_day = fig.add_subplot(gs[-1, 0:2])
# Extract the Pyrano_Combination column for the full day
combinations_full_day = data['Pyrano_Combination']
# Initialize a dictionary to count the occurrences of each combination for the full day
combination_counts_full_day = {}
# Count the occurrences of each combination for the full day
for combo_list in combinations_full_day:
    if isinstance(combo_list, str):
        combo_list = combo_list.replace("'", "").strip("()")
        combos = combo_list.split(", ")
        for combo in combos:
            if combo in combination_counts_full_day:
                combination_counts_full_day[combo] += 1
            else:
                combination_counts_full_day[combo] = 1
# Create a matrix to represent the counts for the heatmap for the full day
heatmap_data_full_day = np.zeros((len(tilt_angles), len(cardinal_directions)))-1
# Populate the heatmap data matrix for the full day
for combo, count in combination_counts_full_day.items():
    if combo != 'GHI':
        direction, angle = combo.split('_')
    else:
        direction, angle = 'GHI', 90
    angle = int(angle)
    if direction in cardinal_directions and angle in tilt_angles:
        dir_index = cardinal_directions.index(direction)
        angle_index = tilt_angles.index(angle)
        heatmap_data_full_day[angle_index, dir_index] = count
# Normalize the heatmap data for the full day
if np.sum(heatmap_data_full_day) > 0:
    heatmap_data_full_day = heatmap_data_full_day / np.sum(heatmap_data_full_day) * 100
heatmap_data_full_day[0,-1]= np.nan
# Plot the heatmap for the full day
c_full_day = ax_full_day.imshow(heatmap_data_full_day, cmap='viridis', aspect='auto', vmin=0, vmax=10)
# Set the labels for the cardinal directions and tilt angles for the full day
ax_full_day.set_xticks(np.arange(len(cardinal_directions)))
ax_full_day.set_xticklabels(cardinal_directions)
ax_full_day.set_yticks(np.arange(len(tilt_angles)))
ax_full_day.set_yticklabels(tilt_angles)
# Set the title for the full day subplot
ax_full_day.set_title('(i) Full Day', fontsize=16, fontweight='bold')
ax_full_day.set_ylabel('Tilt angle [°]')

# Add a single colorbar for the entire figure
cbar_ax = fig.add_subplot(gs[1:4, 2])
cbar = fig.colorbar(c_full_day, cax=cbar_ax)
cbar.set_label(label='%', size=14)  # Adjust the size as needed
cbar.ax.tick_params(labelsize=12)  # Adjust the size as needed

# Adjust layout to prevent overlap
plt.tight_layout(rect=[0, 0, 0.9, 0.95])
# fig.savefig(ALL_PLOT_PATH / "Most_used_pyrano3.png", dpi=300, bbox_inches='tight')
fig.savefig(LOW_RES_PLOT_PATH / "Figure 10.png", dpi=300, bbox_inches='tight')
# fig.savefig(HIGH_RES_PLOT_PATH / "Figure 10.pdf", format='pdf', bbox_inches='tight')
plt.close()
