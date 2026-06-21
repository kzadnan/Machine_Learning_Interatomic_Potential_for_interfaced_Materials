def POSCAR_LAMMPS_ASE(deepmd_order):
    import glob
    import os
    from ase.io import read, write
    from ase.data import atomic_masses, atomic_numbers
    
    # --- 1. DEFINE THE ROBUST CONVERTER FUNCTION ---
    def POSCAR_TO_LAMMPS(input_name, output_name):
        # A. Read the POSCAR
        atoms = read(input_name)
        
        # B. WRITE THE GEOMETRY USING ASE
        write(output_name, atoms, format='lammps-data', 
              specorder=deepmd_order, atom_style='atomic')
        
        # C. DYNAMICALLY GENERATE THE MASS BLOCK
        # This prevents errors by matching the exact deepmd_order provided
        mass_lines = ["\nMasses\n\n"]
        for i, symbol in enumerate(deepmd_order):
            # Fetch the precise mass from ASE's database
            mass = atomic_masses[atomic_numbers[symbol]]
            mass_lines.append(f"{i + 1} {mass:.4f} # {symbol}\n")
        mass_lines.append("\n")
        mass_block = "".join(mass_lines)
        
        # Read the file we just created
        with open(output_name, 'r') as f:
            lines = f.readlines()
            
        # Check if "Masses" exists. If not, insert it.
        if not any("Masses" in line for line in lines):
            insert_idx = 0
            
            # Find insertion point: usually after "xy xz yz" or "zlo zhi"
            for idx, line in enumerate(lines):
                if "xy xz yz" in line:
                    insert_idx = idx + 1
                    break
                    
            if insert_idx == 0:
                for idx, line in enumerate(lines):
                    if "zlo zhi" in line:
                        insert_idx = idx + 1
                        break # FIXED: Indented correctly so it only breaks when found
            
            # Insert the block
            lines.insert(insert_idx, mass_block)
            
            # Write back to file
            with open(output_name, 'w') as f:
                f.writelines(lines)
                
        print(f"converted {input_name} -> {output_name} (Checked Masses & Types)")
    
    # --- 2. MAIN SCRIPT LOOP ---
    file_list = glob.glob("POSCAR-*")
    count = len(file_list)
    
    print(f"Found {count} files. Starting conversion...")
    
    for i in range(1, count + 1):
        input_name = f"POSCAR-{i:03d}"   
        output_name = f"disp-{i:03d}.lammps"
        
        if os.path.exists(input_name):
            POSCAR_TO_LAMMPS(input_name, output_name)
        else:
            print(f"Skipping {input_name} (File not found)")
    
    print("All conversions complete.")