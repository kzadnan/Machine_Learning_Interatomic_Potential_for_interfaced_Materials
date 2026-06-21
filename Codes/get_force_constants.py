# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 16:34:39 2025

@author: Owner
"""

def FORCE_CONSTS(supercell_matrix):
    import numpy as np
    import glob
    from phonopy import load
    from phonopy.file_IO import parse_FORCE_SETS
    from phonopy.file_IO import write_FORCE_CONSTANTS
    from phonopy.interface.vasp import get_vasp_structure_lines
    
    # 1. LOAD THE PHONOPY STATE
    # This re-loads your supercell and displacement info
    print("Loading phonopy_disp.yaml...")
    
    phonon = load("phonopy_disp.yaml",supercell_matrix=supercell_matrix,is_nac=False)
    # Note: modern phonopy.load("phonopy_disp.yaml") is often enough, 
    # but providing the filename ensures safety.
    print(supercell_matrix)
    
    
    
    # 2. DEFINE FUNCTION TO PARSE LAMMPS DUMPS
    def parse_lammps_forces(filename, num_atoms):
        """
        Reads a LAMMPS dump file and returns a (num_atoms, 3) numpy array.
        Assumes columns: id type fx fy fz (or similar, as long as fx,fy,fz are last 3)
        """
        forces = []
        ids = []
        
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Find where the atoms start (Look for "ITEM: ATOMS")
        start_idx = 0
        for i, line in enumerate(lines):
            if "ITEM: ATOMS" in line:
                start_idx = i + 1
                break
                
        # Parse the lines
        for line in lines[start_idx:]:
            data = line.strip().split()
            if len(data) < 4: continue # skip empty lines
            
            # We assume standard dump custom: id type fx fy fz
            # We need the ID to sort them because LAMMPS shuffles atoms!
            atom_id = int(data[0]) 
            fx = float(data[-3])
            fy = float(data[-2])
            fz = float(data[-1])
            
            ids.append(atom_id)
            forces.append([fx, fy, fz])
            
        # Convert to numpy array
        forces = np.array(forces)
        ids = np.array(ids)
        
        # SORT BY ATOM ID
        # This is critical. LAMMPS writes atoms in random order parallelly.
        # Phonopy expects atom 1, atom 2, atom 3... strictly.
        sorted_indices = np.argsort(ids)
        sorted_forces = forces[sorted_indices]
        
        return sorted_forces
    
    # 3. COLLECT FORCES
    # Find all force files. Adjust the pattern to match your filenames!
    # Example: If they are in folders: "disp-*/forces.dump"
    force_files = sorted(glob.glob("disp-*/forces.lammpstrj")) # OR "disp-*/dump.forces"
    
    # If you used the previous scripts, your files might be named "forces.disp-001.lammps"
    if not force_files:
        # Try alternate pattern
        force_files = sorted(glob.glob("forces.disp-*"))
    
    print(f"Found {len(force_files)} force files.")
    
    force_sets = []
    
    # Get number of atoms from the supercell in phonopy
    num_atoms = phonon.supercell.get_number_of_atoms()
    
    for i, filename in enumerate(force_files):
        print(f"Parsing {filename}...")
        f = parse_lammps_forces(filename, num_atoms)
        
        # Check shape
        if f.shape != (num_atoms, 3):
            print(f"❌ Error: {filename} has {len(f)} atoms, expected {num_atoms}!")
            break
            
        force_sets.append(f)
    
    # 4. CALCULATE FORCE CONSTANTS
# ... inside your function ...

# 4. CALCULATE FORCE CONSTANTS
    print("Computing Force Constants...")
    phonon.produce_force_constants(forces=force_sets)
    
    # 5. SAVE TO FILE
    # We use the specific writer from phonopy.file_IO
    from phonopy.file_IO import write_FORCE_CONSTANTS
    
    # Access the calculated constants
    force_constants = phonon.force_constants
    
    # Write to text file
    write_FORCE_CONSTANTS(force_constants, filename="FORCE_CONSTANTS")
    print("✅ Success! FORCE_CONSTANTS generated.")

# ... (Previous code where you saved FORCE_CONSTANTS) ...

    print("Calculating Band Structure...")

    phonon.auto_band_structure(plot=True, write_yaml=True, with_eigenvectors=True)

    # Save band plot
    plt = phonon.plot_band_structure()
    plt.savefig("band_structure.png")

    # Write gnuplot style band.txt
    bands = phonon.get_band_structure_dict()
    
    distances = bands["distances"]      # list of segments
    frequencies = bands["frequencies"]  # list of (nq, nbands)

    nbands = frequencies[0].shape[1]

 
    with open("band.txt", "w") as f:
        for band_index in range(nbands):
            for segment_dist, segment_freq in zip(distances, frequencies):

                for x, freq_vals in zip(segment_dist, segment_freq):
                    f.write(f"{x:.8f} {freq_vals[band_index]:.8f}\n")

                f.write("\n")  # separate segments

            f.write("\n\n")  # separate bands

    print("✅ band.txt written band-by-band.")
    
    # 1. Auto-detect the best path for Diamond/FCC structure
    # (Requires 'seekpath' installed: pip install seekpath)
    # If you don't have seekpath, see Method 2 below.
    
    # 2. Show the plot immediately
   
    plt = phonon.plot_band_structure()
    output_plot = "band_structure.png"
    plt.savefig(output_plot)
    print("✅ Band structure plotted.")


# OPTIONAL: Save HDF5 (Faster for large systems)
# phonon.save(settings={'force_constants': True})
    
if __name__ == "__main__":
    import numpy as np
    supercell_matrix=2*np.diag([1,1,1])
    FORCE_CONSTS(supercell_matrix)
    
    
# OPTIONAL: AUTO-SAVE BAND STRUCTURE CONFIG
# phonon.save("phonopy_params.yaml")