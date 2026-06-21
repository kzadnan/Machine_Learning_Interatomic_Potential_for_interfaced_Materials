# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 12:14:30 2025

@author: Owner
"""
def make_displaced_structures(filename,supercell_matrix,distance):
    import numpy as np
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp, write_vasp
    
    # 1. Load the unit cell
    unitcell = read_vasp("POSCAR")
    
    # 2. Initialize Phonopy with a supercell matrix (e.g., 2x2x2)
    phonon = Phonopy(unitcell,
                     supercell_matrix)
    
    # 3. Generate displacements
    # 'distance=0.01' is the standard displacement distance in Angstroms
    phonon.generate_displacements(distance)
    
    # 4. Get the list of displaced supercells
    # This returns a list of atom objects, one for each required displacement
    displaced_supercells = phonon.supercells_with_displacements
    
    # 5. Write them to files
    for i, d_cell in enumerate(displaced_supercells):
        # We name them POSCAR-001, POSCAR-002, etc. to match Phonopy conventions
        filename = f"POSCAR-{i+1:03d}"
        write_vasp(filename, d_cell, direct=True)
        print(f"Written {filename}")
    
    # Optional: Save 'phonopy_disp.yaml'
    # This file is CRITICAL. It maps the displacement files (POSCAR-001) 
    # to the specific atoms that were moved. You need this later.
    phonon.save("phonopy_disp.yaml")
    print("Written phonopy_disp.yaml")