# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 10:39:11 2025

@author: arthurg
"""
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Define the output directory
output_dir = Path(r"C:/Users/arthurg/OneDrive - Universitetssenteret på Svalbard AS/Documents/UNIS_PhD/PAPER_2/PAPER_2_Data_Analysis/GLOB_scripts/Figures_all")
output_dir.mkdir(exist_ok=True)  # Create directory if it doesn't exist

# Load the data with semicolon separator and set 'Time' as index
file_path = r"C:\Users\arthurg\OneDrive - NTNU\Workspace\Data\Clouds\CloudCoverLYR_2013-2024.csv"
df = pd.read_csv(file_path, sep=";", skiprows=1, parse_dates=["Time"], dayfirst=True, index_col="Time")

# Filter for April to September 2023 and 2024
months = range(4, 10)  # April (4) to September (9)
df_filtered = df[(df.index.month.isin(months)) & ((df.index.year == 2023) | (df.index.year == 2024))]

# Create a new column for 3-hour intervals
df_filtered["Interval"] = (df_filtered.index.hour // 3) * 3

# Create a figure with subplots
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(14, 12))
# fig.suptitle("Difference from Mean Cloud Cover (2023 & 2024) - 3-Hour Intervals", fontsize=16)

# Flatten axes for easy iteration
axes = axes.flatten()

for i, month in enumerate(months):
    # Filter data for the current month
    month_data = df_filtered[df_filtered.index.month == month]
    # Group by the 'Interval' column and calculate the mean cloud cover
    avg_cloud_cover = month_data.groupby("Interval")["Cloud cover"].mean()
    # Calculate the overall mean cloud cover for the month
    overall_mean = avg_cloud_cover.mean()
    # Calculate the difference from the overall mean
    difference_from_mean = avg_cloud_cover - overall_mean

    # Reindex to ensure all intervals are present, even if empty
    all_intervals = pd.Series(index=range(0, 24, 3))
    difference_from_mean = difference_from_mean.reindex(all_intervals.index, fill_value=0)

    # Rotate the data to start at 03:00 and end at 00:00
    rotated_intervals = pd.concat([difference_from_mean.iloc[1:], difference_from_mean.iloc[:1]])

    # Plot bars to the left of the tick
    ax = axes[i]
    bars = ax.bar(rotated_intervals.index, rotated_intervals, width=3, align='edge', color="skyblue")
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(x=12, color='gray', linestyle='--', linewidth=0.8)
    ax.set_title(f"{pd.to_datetime(month, format='%m').month_name()}", fontweight='bold', fontsize=13)
    ax.set_xlabel("Hour of day (CET)")
    ax.set_ylabel("Cloud cover anomaly\n(okta difference)")
    ax.set_ylim(-0.4, 0.4)
    ax.set_xticks(range(0, 25, 3))
    ax.set_xticklabels(['00:00', '03:00', '06:00', '09:00', '12:00', '15:00', '18:00', '21:00', '00:00'])
    ax.set_xlim(0, 24)
    ax.grid(axis="y", linestyle="--", alpha=0.7)
    
    # Remove y-axis labels for all but the rightmost subplots
    if i % 2 != 0:
        ax.set_ylabel("")
    
    # Add annotations only to the first subplot (April)
    if i == 0:
       ax.text(2, 0.35, "more clouds", color="red", ha="center", va="center", fontsize=10)
       ax.text(2, -0.35, "less clouds", color="blue", ha="center", va="center", fontsize=10)

# Adjust layout
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Save the figure
fig.savefig(output_dir / "cloud_cover_difference.png", dpi=300, bbox_inches='tight')
fig.savefig(output_dir/ "Figure high res pdf" / "Figure 11.pdf", bbox_inches='tight')

