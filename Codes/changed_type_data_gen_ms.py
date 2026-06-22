# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 03:25:23 2026

@author: Owner
"""

from dpdata import LabeledSystem, MultiSystems
from glob import glob

# Define your global type map to ensure C-1, N-2, Ga-3, Au-4 mapping
# Note: The order here defines the ID (0, 1, 2, 3)
type_map = ['C', 'N','Al', 'Ti', 'Ga', 'Au']

ms = MultiSystems()
fs = glob("./OUTCAR*")

for f in fs:
    try:
        # Pass the type_map here to force consistent indexing across different OUTCARs
        ls = LabeledSystem(f, fmt='vasp/outcar', type_map=type_map)
        
        if len(ls) > 0:
            ms.append(ls)
            # Correct way to access types in dpdata
            print(f"Successfully added {f} with elements: {ls['atom_names']}")
            
    except Exception as e:
        print(f"Error processing {f}: {e}")


system_names = list(ms.systems.keys())

# Initialize a variable to store the total number of energies
total_energies = 0

# Iterate over each system name
for system_name in system_names:
    # Access the system
    system = ms.systems[system_name]
    
    # Get the number of energy values in this system
    num_energies = len(system.data["energies"])
    
    # Add to the total
    total_energies += num_energies

# Print the total number of energy values
print("Total number of energy values across all systems:", total_energies)

# Save the combined data
# This will create the 'type.raw' file using the IDs from your type_map
#ms.to_deepmd_npy("deepmd_data", type_map=type_map)

#print("\nProcessing complete. 'deepmd_data' is ready for training.")

#%% based on ms.systems

#total=len(ms.systems["C216N0Ga0Au96"].data["energies"]) + len(ms.systems["C216N0Ga0Au96"].data["energies"])
[train_systems, test_systems, test_system_idx]=ms.train_test_split(test_size=int(0.1 * total_energies))

train_systems.to_deepmd_npy('training_data')
test_systems.to_deepmd_npy('validation_data')

#%%

import matplotlib.pyplot as plt

for system_name in system_names:
    plt.plot(ms.systems[system_name].data["energies"],label=system_name)
plt.legend()
