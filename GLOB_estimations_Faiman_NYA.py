# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:35:35 2024

@author: arthurg
Contact: arthurg@unis.no

Description:
This script calculates beam and diffuse irradiance using the Faiman et al. (1992) method.
It processes GLOB and K&Z data, estimates irradiance components, and generates plots.
"""

# %% Load Libraries
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from itertools import combinations
from datetime import datetime

# Add custom script path
sys.path.append(r"C:\Users\arthurg\OneDrive - Universitetssenteret på Svalbard AS\Documents\UNIS_PhD\PAPER_2\PAPER_2_Data_Analysis\GLOB_scripts")

# % Load Data
data_path = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data")

# Load GLOB data
ds_glob = xr.open_dataset(data_path / r"GLOB\GLOB_data_30sec_2025.nc")
lat_glob = ds_glob.latitude.values
lon_glob = ds_glob.longitude.values

# %% Beam and Diffuse Irradiance Calculation
import glob_functions_Faiman as fct

# Define criteria and date range
Criteria = "ERBS"

month = 5  # For example, October
year = 2024

# Create a daily date range for the specified month and year
start_date = f'{year}-{month:02d}-01'
end_date = f'{year}-{month:02d}-{pd.Period(start_date).days_in_month}'
dates = pd.date_range(start=start_date, end=end_date, freq='D')

# Define the azimuth directions and angles
azimuth_directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
angles = [45, 90, 135]
# Variables of interest for irradiance calculation
variables = ['GHI'] + [f"{azimuth}_{angle}" for azimuth in azimuth_directions for angle in angles]


# Generate all combinations of variables
combs = list(combinations(variables, 2))
for date in dates: 
    date=date.strftime('%Y-%m-%d')
    
    df_glob_one_day = ds_glob.sel(Timestamp=date).to_pandas()
    df_glob_one_day.index = df_glob_one_day.index.tz_localize('UTC')

    # Generate all possible combinations
    combs = list(combinations(variables, 2))
        
    # Prepare to store the results
    results = []
    
    # Loop over each minute in the day
    TIMESTAMPS = pd.date_range(start=f"{date} 00:00:00", end=f"{date} 23:50:00", freq='10min', tz='UTC')
    for timestamp in TIMESTAMPS:
        timestamp = pd.DatetimeIndex([timestamp])  # Ensure timestamp is in the desired format
        glob_value = df_glob_one_day.loc[timestamp]
        
        # Calculate the solar position
        solar_position = fct.calculate_solar_angles(timestamp, lat_glob, lon_glob)
        # Extract the zenith angle
        zenith_angle = solar_position['zenith'].values[0]
 
        # Use the optimized function
        D, I, D_prime, I_prime, error, comb_opt = fct.find_best_combination(combs, glob_value, zenith_angle, lat_glob, lon_glob)
        
        comb_opt = str(comb_opt)
        if not np.isnan(D) and not np.isnan(I):
            D = round(D); I = round(I)
            D_prime = round(D_prime); I_prime = round(I_prime)
            error = round(error)
        else: 
            D = str(D); I = str(I)
            D_prime = str(D_prime); I_prime = str(I_prime)
            error = str(np.nan)
        zenith_angle = round(zenith_angle, 2)
        
        # Append the results for this timestamp
        results.append({
            'Timestamp': timestamp[0],
            'Beam': I,  # Beam value
            'Diffuse': D,  # Diffuse value
            'Beam_prime': I_prime,  # Beam value
            'Diffuse_prime': D_prime,  # Diffuse value
            'Error': error, # error = ((I_ref - I)**2 + (D_ref - D)**2)**0.5
            'Best_Combination': comb_opt, # Pyrano combination used for best_estim
            'Solar zenith': zenith_angle,
            'Albedo': round(glob_value['albedo'].values[0], 2)
        })

        print(f"Timestamp: {timestamp[0].strftime('%Y-%m-%d %H:%M:%S')} | I={I} and D={D} W/m2 | Combination:{comb_opt}")

    
    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    output_file = data_path / "GLOB" / "B_and_D_Estimations_NYA" / f"best_estimations_{date}_{Criteria}.csv"
    
    # Write the header and units to the file
    # Get the current date
    current_date = datetime.now().strftime("%Y-%m-%d")
    # Define your name
    your_name = "Arthur Garreau"
    date_production = f"Date of production: {current_date}\n"
    author = f"Produced by: {your_name}\n"
    header = "Best estimation of beam and diffuse irradiance with GLOB using the Faiman et al. (1992) method.\n\
    Location: Ny-Alesund (...N ...E).\n"
    units = "[UTC]\t\t[W m-2]\t [W m-2]\t  [W m-2]\t [W m-2]\t [/]\t [/]\t [°]\t [/]\n"
    
    # Open the file and write the header and units
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(date_production)
        f.write(author)
        f.write(header)
        f.write(units)
    
    # Append the DataFrame to the file
    results_df.to_csv(output_file, index=False, mode='a', sep='\t', encoding='utf-8')


# %% Plot Beam, Diffuse, and Global irradiances estimated by GLOB
import matplotlib.dates as mdates

# dates = pd.date_range(start="2024-05-01", end="2024-05-01", freq="D", tz="UTC")

for date in dates:
        
    date=date.strftime('%Y-%m-%d')
    results_df = pd.read_csv(data_path / "GLOB" / "B_and_D_Estimations_NYA" / f"Faiman_{date}_criteria{Criteria}.csv")  # Adjust path as needed
    results_df['Timestamp'] = pd.to_datetime(results_df['Timestamp'])
    timestamps = results_df['Timestamp']
    
    GHI = ds_glob['GHI'].sel(Timestamp=date).to_pandas()
    zenith = np.deg2rad(results_df['Solar zenith (°)'])
    
    plt.figure(figsize=(9, 6))
    # 
    # Plot LYR Beam
    plt.plot(results_df['Timestamp'], results_df['Beam'], label='Beam from GLOB', color='orange')
    # Plot LYR Diffuse
    plt.plot(results_df['Timestamp'], results_df['Diffuse'], label='Diffuse from GLOB', color='blue')
    # Plot LYR Global
    plt.plot(GHI, label='Global measured with Kipp&Zonen', color='black', alpha=0.5)
    plt.plot(results_df['Timestamp'], np.cos(zenith)*results_df['Beam']+results_df['Diffuse'], 
             label='Global = cos(z).Beam + Diffuse\nfrom GLOB', color='black')
  
    # Customize the plot
    plt.xlabel(f'Time (UTC)\n{date}', fontsize=14); plt.ylabel('Irradiance [W / m²]', fontsize=14)
    plt.legend(fontsize=9)
    plt.grid(True, linestyle=':')
    # plt.tight_layout()

    # Set hourly ticks
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=3))  # Tick every hour
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))  # Format as HH:MM
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    # Show the plot
    save_path = "C:/Users/arthurg/OneDrive - Universitetssenteret på Svalbard AS/Documents/UNIS_PhD/PAPER_2/PAPER_2_Data_Analysis/Fig/Estim_Faiman/"
    plt.savefig(save_path + f"Estim_Faiman_{Criteria}_{date}_NYA")
    # plt.close()
    
# %% Plot Global Vs. reconstructed Global


import numpy as np
import matplotlib.dates as mdates
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

dates = pd.date_range(start="2024-04-14", end="2024-07-30", freq="D", tz="UTC")
df_LYR_ghi = ds_kz["SWdown"].to_dataframe()
df_LYR_ghi = df_LYR_ghi.tz_localize('UTC')
df_LYR_ghi = df_LYR_ghi.sort_index()
# Lists to store data for plotting
all_timestamps = []
all_global_kz = []
all_global_construct = []

# Accumulate data for each date
for date in dates:
    date_str = date.strftime('%Y-%m-%d')
    results_df = pd.read_csv(data_path / "GLOB" / "Best_estimations_Faiman" / f'best_estimations_{date_str}_ERBS.csv')

    # Convert the Timestamp column to datetime
    results_df['Timestamp'] = pd.to_datetime(results_df['Timestamp'])

    # Select the relevant data from df_LYR_ghi
    timestamps = results_df['Timestamp']
    df_LYR_ghi_sel = df_LYR_ghi.loc[timestamps[0]:timestamps.iloc[-1]]
    df_LYR_ghi_sel = df_LYR_ghi_sel.rename(columns={"SWdown": "Global K&Z"})

    # Calculate the irradiance
    zenith_rad = np.deg2rad(results_df['Solar zenith (°)'])
    results_df["Constructed Global"] = np.cos(zenith_rad) * results_df['Beam'] + results_df['Diffuse']

    results_df = pd.merge(df_LYR_ghi_sel, results_df, how='inner', left_on='time', right_on='Timestamp')
    # Append data to lists
    all_global_kz.extend(results_df['Global K&Z'])
    all_global_construct.extend(results_df['Constructed Global'])


# Create a mask for non-NaN values
mask = ~np.isnan(all_global_kz) & ~np.isnan(all_global_construct)
all_global_kz = np.array(all_global_kz); all_global_construct = np.array(all_global_construct)
# Filter the data
all_global_kz = all_global_kz[mask]
all_global_construct = all_global_construct[mask]
# %
# Plot all data
plt.figure(figsize=(12, 8))
# Set font size for all elements
plt.rcParams.update({'font.size': 12})
# Create scatter plot
plt.scatter(all_global_kz, all_global_construct, marker='.', color='k', alpha=0.5)
# Plot the line x = y
plt.plot([min(all_global_kz), max(all_global_kz)], [min(all_global_kz), max(all_global_kz)], color='red', linestyle='--', 
         label='x = y')
# Use Ridge Regression for better numerical stability
model = Ridge(alpha=1e-3)  # alpha is the regularization strength
model.fit(all_global_kz.reshape(-1, 1), all_global_construct)
regression_line = model.predict(all_global_kz.reshape(-1, 1))
coefficients = model.coef_[0]
intercept = model.intercept_
r2 = r2_score(all_global_kz, regression_line)

# Plot the regression line
plt.plot(all_global_kz, regression_line, color='red', 
         label=f'Regression\n y = {coefficients:.2f}.x + {intercept:.2f} | $R^2$: {r2:.2f}')
# Show plot
plt.show()
plt.xlabel('Constructed Global from GLOB (W/m²)')
plt.ylabel('Global from Kipp&Zonen (W/m²)')
plt.title('Data from 2024')
plt.legend()
plt.grid(True, linestyle=':')
plt.show()
save_path = "C:/Users/arthurg/OneDrive - Universitetssenteret på Svalbard AS/Documents/UNIS_PhD/PAPER_2/PAPER_2_Data_Analysis/Fig/"
# plt.savefig(save_path + f"Regression_KZ_GLOB")
# plt.close()