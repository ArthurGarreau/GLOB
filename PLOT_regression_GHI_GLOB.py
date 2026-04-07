# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 17:01:25 2024

@author: arthurg

Contact: arthurg@unis.no
"""
# -*- coding: utf-8 -*-

# %% Data import


import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from config_path import DATA_PATH, FIG_PATH
from sklearn.metrics import mean_squared_error

file_path_glob = DATA_PATH / "GLOB_data_10min_2025.nc"
ds_glob = xr.open_dataset(file_path_glob, engine="h5netcdf")

file_path_kz = DATA_PATH / "NYA_BSRN_data" / "NYA_radiation_2025-all.tab"
df_nya = pd.read_csv(file_path_kz, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time').resample("10min").first()
df_nya.index = df_nya.index.tz_localize('UTC')
df_nya['GHI_nya'] = df_nya['SWD']

file_path_gti = DATA_PATH / "Estim_GTI" / "2025_estimation_GTI_10min_nonlinear_4pyrano.csv"
df_gti = pd.read_csv(file_path_gti, sep='\t', parse_dates=True, index_col='Timestamp', header=10)
df_gti['GHI_estim'] = df_gti['gti0_0']


df_glob = ds_glob.to_dataframe()
df_glob.index = df_glob.index.tz_localize('UTC')
df_glob['GHI_glob'] = df_glob['GHI']

merged_df1 = pd.merge(df_glob['GHI_glob'], df_nya['GHI_nya'], left_index=True, right_index=True, how='inner')
merged_df1 = merged_df1.resample('10min').mean().dropna(); merged_df1 = merged_df1.loc['2025-04-01':]

merged_df2 = pd.merge(df_gti['GHI_estim'], df_nya['GHI_nya'], left_index=True, right_index=True, how='inner').dropna()
merged_df2 = merged_df2.loc['2025-04-01':]

merged_df3 = pd.merge(df_gti['GHI_estim'], df_glob['GHI_glob'], left_index=True, right_index=True, how='inner').dropna()
merged_df3 = merged_df3.loc['2025-04-01':]



# %%

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

# --- Plot for merged_df2 ---
slope2, intercept2, r_value2, p_value2, std_err2 = linregress(merged_df2['GHI_nya'], merged_df2['GHI_estim'])
predicted2 = intercept2 + slope2 * merged_df2['GHI_nya']

rmse2 = np.sqrt(mean_squared_error(merged_df2['GHI_estim'], merged_df2['GHI_nya']))
mbe2 = np.mean(merged_df2['GHI_estim'] - merged_df2['GHI_nya'])

rmse_rounded2 = int(round(rmse2))
mbe_rounded2 = int(round(mbe2))

min_val2 = min(merged_df2['GHI_nya'].min(), merged_df2['GHI_estim'].min())
max_val2 = max(merged_df2['GHI_nya'].max(), merged_df2['GHI_estim'].max())

ax1.plot([min_val2, max_val2], intercept2 + slope2 * np.array([min_val2, max_val2]), 'r', linewidth=2,
         label=f'Linear fit: y={slope2:.2f}x + {intercept2:.2f}\nR² = {r_value2**2:.3f}')
ax1.plot([min_val2, max_val2], [min_val2, max_val2], '--', label='y = x', color='#07BD00', linewidth=2)
ax1.scatter(merged_df2['GHI_nya'], merged_df2['GHI_estim'], label='10-min GHI', color='k', marker='+')

unit = '$W \\ m^{-2}$'
stats_text2 = f'RMSE = {rmse_rounded2} {unit}\nMBE = {mbe_rounded2} {unit}'
ax1.text(0.7, 0.15, stats_text2, transform=ax1.transAxes, verticalalignment='top',
         bbox=dict(facecolor='white', edgecolor='grey', alpha=1, boxstyle='round'))

ax1.set_xlim([-10,900]);ax1.set_ylim([-10,900]);
ax1.set_xlabel(f'GHI Thermopile Ny-Ålesund [{unit}]')
ax1.set_ylabel(f'GHI Estimated [{unit}]')
ax1.legend(edgecolor='none', loc="upper left")
ax1.grid(True, linestyle=":")

# --- Plot for merged_df3 ---
slope3, intercept3, r_value3, p_value3, std_err3 = linregress(merged_df3['GHI_glob'], merged_df3['GHI_estim'])
predicted3 = intercept3 + slope3 * merged_df3['GHI_glob']

rmse3 = np.sqrt(mean_squared_error(merged_df3['GHI_estim'], merged_df3['GHI_glob']))
mbe3 = np.mean(merged_df3['GHI_estim'] - merged_df3['GHI_glob'])

rmse_rounded3 = int(round(rmse3))
mbe_rounded3 = int(round(mbe3))

min_val3 = min(merged_df3['GHI_glob'].min(), merged_df3['GHI_estim'].min())
max_val3 = max(merged_df3['GHI_glob'].max(), merged_df3['GHI_estim'].max())

ax2.plot([min_val3, max_val3], intercept3 + slope3 * np.array([min_val3, max_val3]), 'r', linewidth=2,
         label=f'Linear fit: y={slope3:.2f}x + {intercept3:.2f}\nR² = {r_value3**2:.3f}')
ax2.plot([min_val3, max_val3], [min_val3, max_val3], '--', label='y = x', color='#07BD00', linewidth=2)
ax2.scatter(merged_df3['GHI_glob'], merged_df3['GHI_estim'], label='10-min GHI', color='k', marker='+')

stats_text3 = f'RMSE = {rmse_rounded3} {unit}\nMBE = {mbe_rounded3} {unit}'
ax2.text(0.7, 0.15, stats_text3, transform=ax2.transAxes, verticalalignment='top',
         bbox=dict(facecolor='white', edgecolor='grey', alpha=1, boxstyle='round'))

ax2.set_xlim([-10,900]);ax2.set_ylim([-10,900]);
ax2.set_xlabel(f'GHI Silicon cell GLOB [{unit}]')
ax2.set_ylabel(f'GHI Estimated [{unit}]')
ax2.legend(edgecolor='none', loc="upper left")
ax2.grid(True, linestyle=":")

# Adjust layout and save
plt.tight_layout()
plt.rcParams.update({'font.size': 14})

# plt.savefig(FIG_PATH / "Figure_low_res" / "Figure R1.png", dpi=300, bbox_inches='tight')
# plt.savefig(FIG_PATH / "Figure_high_res" / "Figure R1.pdf", format='pdf', bbox_inches='tight')
# plt.close()

