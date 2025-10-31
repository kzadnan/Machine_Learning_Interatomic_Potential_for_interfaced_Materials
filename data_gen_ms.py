from dpdata import LabeledSystem, MultiSystems
#from glob import glob
"""
process multi systems
"""
import glob

# List all files matching the pattern 'OUTCAR*'
fs = glob.glob("OUTCAR*")

ms = MultiSystems()
for f in fs:
    try:
        ls = LabeledSystem(f,fmt = 'vasp/outcar')
    except:
        print(f)
    if len(ls) > 0:
        ms.append(ls)

#data_all=ms.to_deepmd_raw("deepmd")
#data_all=ms.to_deepmd_npy("deepmd")




#%%
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

#%%

[train_systems, test_systems, test_system_idx]=ms.train_test_split(test_size=int(0.1 * total_energies))

train_systems.to_deepmd_npy('training_data')
test_systems.to_deepmd_npy('validation_data')