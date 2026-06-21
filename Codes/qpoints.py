# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 18:13:36 2025

@author: Khalid
"""
import numpy as np
from phonopy import Phonopy
from phonopy.file_IO import parse_FORCE_CONSTANTS
from phonopy.interface.vasp import read_vasp
import torch

# Load force constants
force_constants_file = "FORCE_CONSTANTS"
force_constants = parse_FORCE_CONSTANTS(force_constants_file)
force_constants_torch = torch.from_numpy(force_constants)


# Load the unit cell structure
unitcell = read_vasp("POSCAR")

# Define the supercell matrix
supercell_matrix = np.eye(3) * 4  # Replace with your actual supercell size

# Initialize Phonopy
phonon = Phonopy(unitcell, supercell_matrix)
phonon.force_constants = force_constants

# Once the dynamical matrix is generated, you can define the q-point mesh
nx=5
ny=nx
nz=nx
mesh = [nx, ny, nz]
phonon.run_mesh(mesh)

phonon.run_mesh(mesh, is_mesh_symmetry=False)

# Extract q-points from the mesh
qpoints = phonon.mesh.qpoints

phonon.run_mesh(mesh, with_group_velocities=True)
#
group_velocity=phonon.group_velocity.group_velocities

dynamical_matrix= phonon.dynamical_matrix.dynamical_matrix

#%% Calculation of dynamical matrix at each qpoint
# Ensure dynamical_matrix_torch has correct dimensions based on your qpoints and system size
num_qpoints = len(qpoints)  # Assuming 4096 here
num_atoms =  len(phonon.unitcell.symbols)# phonon.primitive.get_number_of_atoms()
dynamical_matrix_torch = torch.zeros([num_qpoints, 3 * num_atoms, 3 * num_atoms],dtype=torch.complex128)

# Loop over all q-points to compute the dynamical matrix
for i, qpoint in enumerate(qpoints):
    phonon.dynamical_matrix.run(qpoint)
    dynamical_matrix = phonon.dynamical_matrix.dynamical_matrix
    print(f"Dynamical matrix at q-point {qpoint}:")
    # Ensure torch conversion with correct dtype
    dynamical_matrix_torch[i, :, :] = torch.from_numpy(np.array(dynamical_matrix, dtype=np.complex128))
    
#%% Green's function calculation 

# Parameters
eta = 1e-4  # Convergence parameter
frequency_range = np.linspace(0, 100, 1000)  # Frequency range (in THz or other unit)
identity_matrix = torch.eye(6, dtype=torch.complex128)

# Initialize storage for Green's functions
greens_functions_torch = torch.zeros(
    (len(qpoints), len(frequency_range), 3 * num_atoms, 3 * num_atoms), dtype=torch.complex128
)

#%%
# Compute Green's functions
# for i, dynamical_matrix in enumerate(dynamical_matrix_torch):
#     for j, omega in enumerate(frequency_range):
#         omega_complex = torch.tensor(omega + 1j * eta, dtype=torch.complex128)  # Add small imaginary component
#         omega_squared_I = (omega_complex**2) * identity_matrix

#         # Compute the Green's function
#         try:
#             greens_functions_torch[i, j, :, :] = torch.linalg.inv(
#                 omega_squared_I - dynamical_matrix
#             )
#         except torch.linalg.LinAlgError:
#             print(f"Matrix inversion failed at q-point {i} and frequency {omega} THz.")


#%%

# import matplotlib.pyplot as plt

# q_idx = 34  # Select the q-point index
# freq_idx = 10# Select the frequency index

# # Assuming the Green's function is 6x6 for each q-point and frequency
# gf_matrix = greens_functions_torch[q_idx, freq_idx, :, :].numpy()

# # Plot the real and imaginary parts of the Green's function
# plt.figure(figsize=(20, 5))

# plt.subplot(1, 2, 1)
# plt.imshow(np.real(gf_matrix), cmap='viridis')
# plt.title(f'Real part of Green\'s\n function at q-point {qpoints[q_idx]} and frequency {frequency_range[freq_idx]}')
# plt.colorbar()

# plt.subplot(1, 2, 2)
# plt.imshow(np.imag(gf_matrix), cmap='viridis')
# plt.title(f'Imaginary part of Green\'s function at q-point {qpoints[q_idx]} and frequency {frequency_range[freq_idx]}')
# plt.colorbar()

# plt.show()



