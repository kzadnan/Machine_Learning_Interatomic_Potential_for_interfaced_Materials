# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 11:57:40 2025

@author: Khalid Zobaid Adnan
"""



import numpy as np
from make_sposcar import make_sposcar
from make_displaced_structures import make_displaced_structures
from Lammps_poscar import Lammps_poscar
from poscar_lammps_ase import POSCAR_LAMMPS_ASE
from organize import organize_folders
from lmp_run import running_lammps
from get_force_constants import FORCE_CONSTS
from get_properties_gpt import properties

#Concerting the Conventional cell to POSCAR format
species_map=  {1: 'Au'}

deepmd_order = [species_map[i] for i in sorted(species_map.keys())]

 #C480N0Al0Ti0Ga0Au0
Lammps_poscar("Au_final_relaxed_structure.txt","POSCAR",species_map)
## Defining the supercell size 2*2*2
supercell_matrix=np.array([[4,0,0],[0,4,0],[0,0,4]])#4*np.diag([1,1,1])
## Making SPOSCAR
make_sposcar("POSCAR", supercell_matrix)
## Making POSCAR-001 -- From irreducible minimum perturbations
displacement=0.2
make_displaced_structures("POSCAR",supercell_matrix,displacement)
##Converting the POSCARs to LAMMPS format
POSCAR_LAMMPS_ASE(deepmd_order)
## Organize files for Force calculations
organize_folders()
## running lammps with deepmd is the issue for each folder in Windows.
running_lammps()
##Getting the force constants
FORCE_CONSTS(supercell_matrix)
#%% Getting running kmesh calculations for phonopy 
mesh=[16,16,16]
properties(supercell_matrix, mesh)







