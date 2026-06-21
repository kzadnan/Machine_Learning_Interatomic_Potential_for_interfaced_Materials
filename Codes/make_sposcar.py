# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 11:54:19 2025

@author: Owner
"""

#%%create SPOSCAR from POSCAR

def make_sposcar(filename,supercell_matrix):
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp, write_vasp
    
    # 1. Load the unit cell from POSCAR
    unitcell = read_vasp(filename)
    
    # 2. Define the Supercell Matrix
    # Example: A 2x2x2 supercell
    # supercell_matrix = [[2, 0, 0],
    #                     [0, 2, 0],
    #                     [0, 0, 2]]
    
    # 3. Instantiate the Phonopy object
    # This automatically generates the supercell internally based on the matrix
    phonon = Phonopy(unitcell,
                     supercell_matrix=supercell_matrix)
    
    # 4. Get the supercell structure object
    supercell = phonon.supercell
    
    # 5. Write the supercell to a file named 'SPOSCAR'
    # direct=True ensures coordinates are in direct (fractional) coordinates
    write_vasp("SPOSCAR", supercell, direct=True)
    print("SPOSCAR generated successfully.")
    return supercell
