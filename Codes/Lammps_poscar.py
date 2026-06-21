# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 12:48:14 2025

@author: Owner
"""
def Lammps_poscar(input_lammps,output_poscar,species_map):
    import numpy as np
    from ase.io import read, write
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp, write_vasp
    
    # ==========================================
    # CONVERT LAMMPS -> POSCAR
    # ==========================================
    
    # input_lammps = "Cunitcell.txt"
    # output_poscar = "POSCAR"
    
    print(f"Reading {input_lammps}...")
    
    # 1. Read the LAMMPS file
    # We use style='atomic' because your file format matches "id type x y z"
    atoms = read(input_lammps, format='lammps-data', atom_style='atomic')
    
    # 2. FORCE CORRECT CHEMICAL SYMBOLS
    # Your file says: Type 1 = Si, Type 2 = C.
    # We must map the integer types (1, 2) to chemical symbols ('Si', 'C').
    
    # Get the list of types from the file (returns array like [2, 2, 1, ...])
    # Note: ASE reads LAMMPS types into the 'type' array.
    types = atoms.get_array('type')
    
    # Create a list of symbols based on your definitions
    # Type 1 (index 0) -> Si
    # Type 2 (index 1) -> C
    #species_map = {1: 'Ga', 2:'N', 3: 'C'}
    
    # Generate the symbol list for every atom in the file
    chemical_symbols = [species_map[t] for t in types]
    
    # Apply these symbols to the atoms object
    atoms.set_chemical_symbols(chemical_symbols)
    
    # 3. Write the POSCAR file
    write(output_poscar, atoms, format='vasp', direct=True, vasp5=True)
    print(f"✅ Converted to {output_poscar} (mapped Type 1->chemical_symbols[0], Type 2->chemical_symbols[1])")