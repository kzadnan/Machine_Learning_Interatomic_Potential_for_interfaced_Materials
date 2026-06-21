# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 03:46:59 2025

@author: Owner
"""

import yaml
import numpy as np

filename = "band.yaml"

print(f"Reading {filename}... (This may take a moment for large files)")
try:
    # use CLoader for speed if available, otherwise safe_load
    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader

    with open(filename, 'r') as f:
        data = yaml.load(f, Loader=Loader)

    # 1. Extract Frequencies (Eigenvalues)
    # The structure is: data['phonon'] -> list of q-points -> 'band' -> list of bands -> 'frequency'
    all_frequencies = []
    distances = []
    
    for q_point in data['phonon']:
        dist = q_point['distance']
        # Get list of frequencies for this specific q-point
        freqs = [band['frequency'] for band in q_point['band']]
        
        distances.append(dist)
        all_frequencies.append(freqs)

    # Convert to numpy array
    # Shape will be (Number_of_qpoints, Number_of_bands)
    eigenvalues = np.array(all_frequencies)
    distances = np.array(distances)

    print(f"✅ Extracted eigenvalues for {eigenvalues.shape[0]} q-points and {eigenvalues.shape[1]} bands.")

    # 2. (Optional) Save to a text file
    # Col 1: Distance, Cols 2+: Eigenvalues for each band
    output_data = np.column_stack((distances, eigenvalues))
    np.savetxt("eigenvalues.txt", output_data, fmt='%15.8f', 
               header="Distance  Band_1  Band_2  Band_3 ...")
    print("✅ Saved to eigenvalues.txt")

except FileNotFoundError:
    print("❌ Error: band.yaml not found.")