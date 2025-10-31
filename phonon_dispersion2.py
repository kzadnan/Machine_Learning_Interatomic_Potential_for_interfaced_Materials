# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:29:50 2025

@author: Khalid
"""

import numpy as np
import matplotlib.pyplot as plt

# Function to load data from a band.txt file


def load_band_data(file_path):
    distances = []
    frequencies = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines or comment lines starting with '#'
            if not line or line.startswith('#'):
                continue

            # Parse the numeric data
            data = list(map(float, line.split()))
            distances.append(data[0])  # First column is the wave vector distance
            frequencies.append(data[1:])  # Remaining columns are frequencies

    distances = np.array(distances)
    frequencies = np.array(frequencies)
    return distances, frequencies


# Load data from the two band.txt files
distances1, frequencies1 = load_band_data("band_dft.txt")
distances2, frequencies2 = load_band_data("band_bulk_mtp.txt")
distances3, frequencies3 = load_band_data("band_interface_deepmd.txt")
distances4, frequencies4 = load_band_data("band_interface_mtp.txt")
distances5, frequencies5 = load_band_data("band_Tersoff.txt")

# Plot settings
plt.figure(figsize=(8, 6))

# Plot the first dataset
for i in range(frequencies1.shape[1]):
    plt.plot(distances1, frequencies1[:, i], '*', label='DFT-Diamond' if i == 0 else "")

# Plot the second dataset
for i in range(frequencies2.shape[1]):
    plt.plot(distances2, frequencies2[:, i], 'r*', label='MTP-Dia' if i == 0 else "")


for i in range(frequencies3.shape[1]):
    plt.plot(distances3, frequencies3[:, i], 'g*', label='DeepMD-Si-Dia' if i == 0 else "")

for i in range(frequencies4.shape[1]):
    plt.plot(distances4, frequencies4[:, i], 'c*', label='MTP-Si-Dia' if i == 0 else "")
    
for i in range(frequencies5.shape[1]):
    plt.plot(distances5, frequencies5[:, i], 'm*', label='Tersoff-dia' if i == 0 else "")    
    
# Add labels, legend, and title
plt.xlabel("Wave Vector")
plt.ylabel("Frequency (THz)")
plt.title("Phonon Dispersion Curve of Diamond")
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save the plot or show it
plt.savefig("phonon_dispersion_comparison.png", dpi=300)
plt.show()
