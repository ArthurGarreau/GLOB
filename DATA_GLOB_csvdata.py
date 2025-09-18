# -*- coding: utf-8 -*-
"""
GLOB Data Processing Script
===========================

This script processes GLOB data from different years (2023, 2024, 2025) and locations (LYR, NYA).
It performs the following tasks:
1. Loads raw GLOB data from CSV files.
2. Renames columns to standardize naming conventions.
3. Converts timestamps to UTC.
4. Resamples data to 30-second intervals.
5. Replaces negative values with a placeholder (-9999).
6. Writes the processed data to new CSV files, preserving specific header lines from the original files.

Key Features:
-------------
- Handles timezone conversions for different periods.
- Ensures data consistency by standardizing column names and time intervals.
- Preserves important metadata from the original data files.

Dependencies:
-------------
- pandas
- pathlib

Author: Arthur Garreau
Contact: arthurg@unis.no
Date: November 1, 2024
"""

import pandas as pd
from pathlib import Path


############################## File Paths #####################################

RAW_DATA = Path(r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\GLOB\Raw_Data")
input_file_data_path_2023 = RAW_DATA / "DATA_2023"
input_file_data_path_2024 = RAW_DATA / "DATA_2024"
input_file_data_path_2025 = RAW_DATA / "DATA_2025"

output_file_data_path = RAW_DATA.parent
###############################################################################

"""
# %% Data 30 seconds 2023
# Load the data into a DataFrame
input_file  = input_file_data_path_2023 / "CR1000X GLOB_Sec.dat" 
output_file = output_file_data_path / "GLOB_data_30sec_2023_LYR.dat"

# Load the data, skipping the header rows for processing
df = pd.read_csv(input_file, skiprows=[0,2,3])

# Rename columns as specified
df = df.rename(columns={
    'SIR(1)': 'GHI',
    'SIR(2)': 'S_45',
    'SIR(3)': 'S_90',
    'SIR(4)': 'S_135',
    'SIR(5)': 'SW_45',
    'SIR(6)': 'SW_90',
    'SIR(7)': 'SW_135',
    'SIR(8)': 'W_45',
    'SIR(9)': 'W_90',
    'SIR(10)': 'W_135',
    'SIR(11)': 'NW_45',
    'SIR(12)': 'NW_90',
    'SIR(13)': 'NW_135',
    'SIR(14)': 'N_45',
    'SIR(15)': 'N_90',
    'SIR(16)': 'N_135',
    'SIR(17)': 'NE_45',
    'SIR(18)': 'NE_90',
    'SIR(19)': 'NE_135',
    'SIR(20)': 'E_45',
    'SIR(21)': 'E_90',
    'SIR(22)': 'E_135',
    'SIR(23)': 'SE_45',
    'SIR(24)': 'SE_90',
    'SIR25': 'SE_135',
    'Kglob': 'GHI_ground',
    'TIMESTAMP': 'Timestamp'
})


# Convert Timestamp to datetime and set it as index
df['Timestamp'] = pd.to_datetime(df['Timestamp'])
df = df.set_index('Timestamp')

# Ensure the index is a DatetimeIndex
df.index = pd.DatetimeIndex(df.index)

# Split data based on the specified date for timezone localization
before_split = df[df.index < '2023-06-20']
after_split = df[df.index >= '2023-06-20']

before_split.index = before_split.index.tz_localize('Etc/GMT-1').tz_convert('UTC')
after_split.index = after_split.index.tz_localize('Etc/GMT-4').tz_convert('UTC')

# Concatenate both parts back together
df = pd.concat([before_split, after_split])

# Resample the data to every 30 seconds, keeping only one sample per interval
df_30s = df.resample('30S').first()

Temp = df_30s['Temp']

df_30s = df_30s.where(df_30s >= 0, -9999)
df_30s = df_30s.fillna(-9999).round().astype(int)
df_30s = df_30s.round().astype(int)

df_30s['Temp'] = Temp.round(1)

df_30s = df_30s[df_30s.index > '2023-02-23 12:00:00+00:00']


df_30s.to_csv(output_file, header=True, index=True)


# Read specific lines from the input file
with open(input_file, 'r') as infile:
    lines = infile.readlines()
    selected_lines = [lines[i] for i in [0, 2, 3]]  # Lines 0, 2, and 3

# Write selected lines at the beginning of the output file
with open(output_file, 'r+') as outfile:
    # Read the existing content
    existing_content = outfile.read()
    
    # Move cursor to the beginning of the file to write new lines
    outfile.seek(0)
    
    # Write the selected lines followed by the existing content
    outfile.writelines(selected_lines)
    outfile.write(existing_content)

df= pd.read_csv(output_file)


# %% Data 30 seconds 2024

import pandas as pd

# Load the data into a DataFrame
input_file  = input_file_data_path_2024 / "CR1000X GLOB_Sec.dat" 
output_file = output_file_data_path / "GLOB_data_30sec_2024_LYR.dat"


# Load the data, skipping the header rows for processing
df = pd.read_csv(input_file, skiprows=[0,2,3], low_memory=False)

# Rename columns as specified
df = df.rename(columns={
    'SIR(1)': 'GHI',
    'SIR(2)': 'S_45',
    'SIR(3)': 'S_90',
    'SIR(4)': 'S_135',
    'SIR(5)': 'SW_45',
    'SIR(6)': 'SW_90',
    'SIR(7)': 'SW_135',
    'SIR(8)': 'W_45',
    'SIR(9)': 'W_90',
    'SIR(10)': 'W_135',
    'SIR(11)': 'NW_45',
    'SIR(12)': 'NW_90',
    'SIR(13)': 'NW_135',
    'SIR(14)': 'N_45',
    'SIR(15)': 'N_90',
    'SIR(16)': 'N_135',
    'SIR(17)': 'NE_45',
    'SIR(18)': 'NE_90',
    'SIR(19)': 'NE_135',
    'SIR(20)': 'E_45',
    'SIR(21)': 'E_90',
    'SIR(22)': 'E_135',
    'SIR(23)': 'SE_45',
    'SIR(24)': 'SE_90',
    'SIR25': 'SE_135',
    'Kglob': 'GHI_ground',
    'TIMESTAMP': 'Timestamp'
})


# Convert Timestamp to datetime and set it as index
df['Timestamp'] = pd.to_datetime(df['Timestamp'])
df = df.set_index('Timestamp')

# Ensure the index is a DatetimeIndex
df.index = pd.DatetimeIndex(df.index)


df.index = df.index.tz_localize('Etc/GMT-4').tz_convert('UTC')

df_30s = df

Temp = df_30s['Temp']

df_30s = df_30s.where(df_30s >= 0, -9999)
df_30s = df_30s.fillna(-9999).round().astype(int)
df_30s = df_30s.round().astype(int)

df_30s['Temp'] = Temp.round(1)
df_30s = df_30s[df_30s.index > '2023-02-23 12:00:00+00:00']


df_30s.to_csv(output_file, header=True, index=True)


# Read specific lines from the input file
with open(input_file, 'r') as infile:
    lines = infile.readlines()
    selected_lines = [lines[i] for i in [0, 2, 3]]  # Lines 0, 2, and 3

# Write selected lines at the beginning of the output file
with open(output_file, 'r+') as outfile:
    # Read the existing content
    existing_content = outfile.read()
    
    # Move cursor to the beginning of the file to write new lines
    outfile.seek(0)
    
    # Write the selected lines followed by the existing content
    outfile.writelines(selected_lines)
    outfile.write(existing_content)

df= pd.read_csv(output_file, low_memory=False)
"""
# %% Data 30 seconds 2025 NYA

# Load the data into a DataFrame
input_file  = input_file_data_path_2025 / "CR1000X_GLOB_online_GLOB_30Sec.dat" 
output_file = output_file_data_path / "GLOB_data_30sec_2025_NYA.dat"


# Load the data, skipping the header rows for processing
df = pd.read_csv(input_file, skiprows=[0,2,3], low_memory=False)

# Rename columns as specified
df = df.rename(columns={
    'SIR(1)': 'GHI',
    'SIR(2)': 'S_45',
    'SIR(3)': 'S_90',
    'SIR(4)': 'S_135',
    'SIR(5)': 'SW_45',
    'SIR(6)': 'SW_90',
    'SIR(7)': 'SW_135',
    'SIR(8)': 'W_45',
    'SIR(9)': 'W_90',
    'SIR(10)': 'W_135',
    'SIR(11)': 'NW_45',
    'SIR(12)': 'NW_90',
    'SIR(13)': 'NW_135',
    'SIR(14)': 'N_45',
    'SIR(15)': 'N_90',
    'SIR(16)': 'N_135',
    'SIR(17)': 'NE_45',
    'SIR(18)': 'NE_90',
    'SIR(19)': 'NE_135',
    'SIR(20)': 'E_45',
    'SIR(21)': 'E_90',
    'SIR(22)': 'E_135',
    'SIR(23)': 'SE_45',
    'SIR(24)': 'SE_90',
    'SIR25': 'SE_135',
    'Kglob': 'GHI_ground',
    'TIMESTAMP': 'Timestamp'
})


# Convert Timestamp to datetime and set it as index
df['Timestamp'] = pd.to_datetime(df['Timestamp'])
df = df.set_index('Timestamp')

# Ensure the index is a DatetimeIndex
df.index = pd.DatetimeIndex(df.index)


df.index = df.index.tz_localize('UTC')

df_30s = df

Temp = df_30s['Temp']

df_30s = df_30s.where(df_30s >= 0, -9999)
df_30s = df_30s.fillna(-9999).round().astype(int)
df_30s = df_30s.round().astype(int)

df_30s['Temp'] = Temp.round(1)


df_30s.to_csv(output_file, header=True, index=True)


# Read specific lines from the input file
with open(input_file, 'r') as infile:
    lines = infile.readlines()
    selected_lines = [lines[i] for i in [0, 2, 3]]  # Lines 0, 2, and 3

# Write selected lines at the beginning of the output file
with open(output_file, 'r+') as outfile:
    # Read the existing content
    existing_content = outfile.read()
    
    # Move cursor to the beginning of the file to write new lines
    outfile.seek(0)
    
    # Write the selected lines followed by the existing content
    outfile.writelines(selected_lines)
    outfile.write(existing_content)

df= pd.read_csv(output_file, low_memory=False)






