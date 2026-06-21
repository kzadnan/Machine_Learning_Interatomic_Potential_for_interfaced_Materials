def properties(supercell_matrix, mesh):
    import numpy as np
    import os
    from phonopy import Phonopy
    from phonopy.file_IO import parse_FORCE_CONSTANTS
    from phonopy.interface.vasp import read_vasp
    
    # Load force constants
    force_constants_file = "FORCE_CONSTANTS"
    force_constants = parse_FORCE_CONSTANTS(force_constants_file)
    
    # Load the unit cell structure
    unitcell = read_vasp("POSCAR")
    
    # Initialize Phonopy
    phonon = Phonopy(unitcell, supercell_matrix)
    phonon.force_constants = force_constants
    
    # Extract properties from the mesh
    phonon.run_mesh(mesh, with_group_velocities=True, with_eigenvectors=True, is_mesh_symmetry=False)
    phonon.write_yaml_mesh()
    
    mesh_dict = phonon.get_mesh_dict()
    qpoints = mesh_dict['qpoints']
        
    np.savetxt('QPOINTS', qpoints, fmt='%.6f')
    np.savetxt('eigenvalues.txt', mesh_dict['frequencies'], fmt='%.6f')
    
    # --- CORRECTED EIGENVECTOR EXTRACTION ---
    eigenvec = mesh_dict['eigenvectors']
    flattened_eigenvec = np.real(eigenvec).reshape(eigenvec.shape[0], -1)
    np.savetxt('eigenvectors_real.txt', flattened_eigenvec, fmt='%.6f')
    print("✅ Properties successfully extracted and saved.")
    
    # =========================================================
    # --- ANIMATE HIGH-SYMMETRY POINTS (UNIT CELL + LATTICE) ---
    # =========================================================
    
    target_qpoints = {
        "Gamma": [0.0, 0.0, 0.0],
        "X":     [0.5, 0.0, 0.5],
        "L":     [0.5, 0.5, 0.5]
    }
    
    num_bands = mesh_dict['frequencies'].shape[1]
    
    # CRITICAL FIX: Extract the lattice matrix of the UNIT CELL (4 atoms), 
    # not the supercell, so the box perfectly wraps the exported atoms.
    cell_matrix = phonon.unitcell.get_cell()
    lattice_str = 'Lattice="' + ' '.join([f"{x:.6f}" for x in cell_matrix.flatten()]) + '"'
    
    for point_name, current_q_point in target_qpoints.items():
        print(f"Generating animations for the {point_name} point...")
        
        folder_name = f"animations_{point_name}"
        os.makedirs(folder_name, exist_ok=True)
        
        for band in range(1, num_bands + 1):
            filename = os.path.join(folder_name, f"{point_name}_band_{band:04d}.xyz")
            
            # REMOVED the invalid 'cell' argument. It now safely outputs the 4-atom cell.
            phonon.write_animation(
                q_point=current_q_point, 
                anime_type='xyz', 
                band_index=band,       
                amplitude=5.0, 
                num_div=20, 
                filename=filename
            )
            
            # Inject the unit cell Lattice Dimensions for OVITO
            with open(filename, 'r') as f:
                lines = f.readlines()
                
            if not lines: continue
                
            num_atoms = int(lines[0].strip())
            step_size = num_atoms + 2
            
            for i in range(1, len(lines), step_size):
                original_comment = lines[i].strip()
                lines[i] = f"{original_comment} {lattice_str}\n"
                
            with open(filename, 'w') as f:
                f.writelines(lines)
                
    print("✅ All high-symmetry animations saved safely with embedded unit box dimensions!")