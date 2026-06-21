import numpy as np
from phonopy import Phonopy
from phonopy.file_IO import parse_FORCE_CONSTANTS
from phonopy.interface.vasp import read_vasp

# Load force constants
force_constants_file = "FORCE_CONSTANTS"
force_constants = parse_FORCE_CONSTANTS(force_constants_file)

# Load the unit cell structure
unitcell = read_vasp("POSCAR")

# Define the supercell matrix
supercell_matrix = np.eye(3) * 2  # Replace with your actual supercell size

# Initialize Phonopy
phonon = Phonopy(unitcell, supercell_matrix)
phonon.set_force_constants(force_constants)

# Define the list of q-points (e.g., Gamma and some neighboring q-points)
q_points = [
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0.5, 0.0],  # Some other q-point
    [0.5, 0.0, 0.0],  # Another q-point
    [0.0, 0.5, 0.5]   # Another q-point
]

# Run phonon calculations for the given q-points
phonon.run_qpoints(q_points=q_points, with_eigenvectors=True)

# Extract phonon data for all q-points
qpoints_dict = phonon.get_qpoints_dict()

# Loop through each q-point and print frequencies and eigenvectors
for idx, q_point in enumerate(q_points):
    frequencies = qpoints_dict['frequencies'][idx]
    eigenvectors = qpoints_dict['eigenvectors'][idx]

    print(f"Q-point {q_point}:")
    print("Frequencies (THz):", frequencies)
    print("Eigenvectors:")
    print(eigenvectors)
    print("-" * 30)
