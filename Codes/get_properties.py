# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 23:57:49 2025

@author: Owner
"""

def properties(supercell_matrix,mesh):
    import numpy as np
    from phonopy import Phonopy
    from phonopy.file_IO import parse_FORCE_CONSTANTS
    from phonopy.interface.vasp import read_vasp
    #import torch
    
    # Load force constants
    force_constants_file = "FORCE_CONSTANTS"
    force_constants = parse_FORCE_CONSTANTS(force_constants_file)
    #force_constants_torch = torch.from_numpy(force_constants)
    
    
    # Load the unit cell structure
    unitcell = read_vasp("POSCAR")
    
    # Define the supercell matrix
    #supercell_matrix = np.eye(3) * 2  # Replace with your actual supercell size
    
    # Initialize Phonopy
    phonon = Phonopy(unitcell, supercell_matrix)
    phonon.force_constants = force_constants
    

    # Extract properties from the mesh
    
    phonon.run_mesh(mesh, with_group_velocities=True,with_eigenvectors=True,
                    
                    is_mesh_symmetry=False)
    
    phonon.write_yaml_mesh()
    #
    #group_velocity=phonon.group_velocity.group_velocities
    
    #dynamical_matrix= phonon.dynamical_matrix.dynamical_matrix
    
    mesh_dict = phonon.get_mesh_dict()
    
    
    qpoints = mesh_dict['qpoints']
    
        
    np.savetxt('QPOINTS', qpoints,fmt='%.6f')
    
    np.savetxt('eigenvalues.txt', mesh_dict['frequencies'], fmt='%.6f')
    eigenvec=mesh_dict['eigenvectors'][0]
    for i in range(1,mesh_dict['eigenvectors'].shape[0]):
        eigenvec=np.vstack([eigenvec, mesh_dict['eigenvectors'][i]])
    split_data = eigenvec.view(float)
    np.savetxt('eigenvectors.txt', split_data, fmt='%.6f')

