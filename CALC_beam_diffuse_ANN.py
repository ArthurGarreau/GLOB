# -*- coding: utf-8 -*-
"""
Created on Mon May 19 16:29:31 2025

@author: arthurg
"""

from config_path import SCRIPT_PATH, DATA_PATH
import pandas as pd
import xarray as xr
import numpy as np
import tensorflow
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential, load_model
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# %% functions definition

def load_and_preprocess_bsrn_data(file_path, frequency):
    """
    Load and preprocess BSRN data from a file.
    """
    bsrn_data = pd.read_csv(file_path, sep='\t', skiprows=24, parse_dates=['Date/Time'], index_col='Date/Time')
    bsrn_data = bsrn_data.resample(f'{frequency}min').first()
    bsrn_data = bsrn_data[['DIF [W/m**2]', 'DIR [W/m**2]']]
    return bsrn_data

def load_and_preprocess_glob_data(file_path):
    """
    Load and preprocess GLOB data from a NetCDF file.
    """
    ds_glob = xr.open_dataset(file_path)
    ds_df_glob = ds_glob.to_dataframe()
    return ds_df_glob

def merge_datasets(bsrn_data, glob_data):
    """
    Merge BSRN and GLOB datasets on the index (timestamps).
    """
    merged_df = pd.merge(glob_data, bsrn_data, left_index=True, right_index=True, how='inner')
    merged_df = merged_df.dropna()
    return merged_df

def prepare_features_and_targets(merged_df):
    """
    Prepare input features and target variables.
    """
    X = merged_df.drop(columns=['DIF [W/m**2]', 'DIR [W/m**2]', 'RECORD','Temp']).values
    y = merged_df[['DIF [W/m**2]', 'DIR [W/m**2]']].values
    return X, y

def split_and_standardize_data(X, y):
    """
    Split the data into training and validation sets and standardize the input features.
    """
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    return X_train, X_val, y_train, y_val, scaler

def define_and_compile_model(input_shape):
    """
    Define and compile the neural network model.
    """
    model = keras.Sequential([
        layers.Dense(128, activation='relu', input_shape=input_shape),
        layers.Dense(64, activation='relu'),
        layers.Dense(2)  # Output layer for DIR and DIF
    ])
    model.compile(optimizer='adam', loss='mse', metrics=['mae'])
    return model

def train_and_evaluate_model(model, X_train, y_train, X_val, y_val, epochs=10, batch_size=32):
    """
    Train and evaluate the neural network model.
    """
    history = model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, validation_data=(X_val, y_val))
    val_loss, val_mae = model.evaluate(X_val, y_val)
    print(f"Validation Loss: {val_loss}")
    print(f"Validation MAE: {val_mae}")
    return history

# %% Main code for the ANN training 


f = 1  # min data frequency
year = 2025

bsrn_datafile = DATA_PATH.parent / "Irradiance" / "NYA" / "NYA_radiation_2025-03.tab"
glob_nc_datafile = DATA_PATH / f"GLOB_data_{f}min_{year}.nc"

bsrn_data = load_and_preprocess_bsrn_data(bsrn_datafile, f)
glob_data = load_and_preprocess_glob_data(glob_nc_datafile)
glob_data = glob_data.reindex(bsrn_data.index)

merged_df = merge_datasets(bsrn_data, glob_data)
X, y = prepare_features_and_targets(merged_df)
X_train, X_val, y_train, y_val, scaler = split_and_standardize_data(X, y)

model = define_and_compile_model((X_train.shape[1],))
history = train_and_evaluate_model(model, X_train, y_train, X_val, y_val)

# Save the model
model.save(SCRIPT_PATH / 'Data' / 'model_estimation_BandI_ANN.keras')


# %% Estimate and Plot

import pvlib
import matplotlib.pyplot as plt

loaded_model = load_model(SCRIPT_PATH / 'Data' / 'model_estimation_BandI_ANN.keras')

# Define the date range for prediction
start_date = '2025-03-29'
end_date = '2025-04-05'
date_range = pd.date_range(start=start_date, end=end_date, freq=f'{f}min')

# Filter the ds data for the specified date range
df_glob_sel = pd.DataFrame(glob_data, index=date_range).dropna()

# Drop rows with NaN values
df_glob_sel = df_glob_sel.drop(columns=['RECORD','Temp'])

# Standardize the April input features using the same scaler
X_april = scaler.transform(df_glob_sel)

# Make predictions using the trained model
predictions = loaded_model.predict(X_april)
diffuse, direct = predictions[:, 0], predictions[:, 1]

# Calculate Diffuse and Direct Irradiance using ERBS model
# Assuming you have the necessary parameters for the ERBS model
# For example, you need total irradiance (ghi) and other parameters
ghi = df_glob_sel['GHI'].values  # Assuming 'GHI' is the column name for global horizontal irradiance
solar_zenith = df_glob_sel['zenith'].values  # Assuming 'zenith_angle' is the column name for solar zenith angle

# Use the ERBS model to split GHI into DNI and DHI
ERBS = pvlib.irradiance.erbs(ghi, solar_zenith, df_glob_sel.index)
dni_erbs, dhi_erbs = ERBS['dni'], ERBS['dhi'] 

bsrn_data_sel = bsrn_data[start_date:end_date]

# Plot the predictions and ERBS model results
plt.figure(figsize=(12, 6))

# Plot Diffuse Irradiance
plt.plot(df_glob_sel.index, diffuse, "r:", label='Neural Network Diffuse Irradiance')
plt.plot(df_glob_sel.index, dhi_erbs,"b:", label='ERBS Diffuse Irradiance')
plt.plot(bsrn_data_sel.index, bsrn_data_sel["DIF [W/m**2]"], "k:", label='BSRN NyÅ Diffuse')

# Plot Direct Irradiance
plt.plot(df_glob_sel.index, direct, "r", label='Neural Network Direct Irradiance')
plt.plot(df_glob_sel.index, dni_erbs, "b" ,label='ERBS Direct Irradiance')
plt.plot(bsrn_data_sel.index, bsrn_data_sel["DIR [W/m**2]"], "k", label='BSRN NyÅ Direct')

# Add labels and title
plt.xlabel('Date')
plt.ylabel('Irradiance (W/m²)')
plt.legend()
# Set hourly ticks
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim((0,800))
plt.grid(True, linestyle=':')

# Show the plot
plt.grid(True)
plt.show()