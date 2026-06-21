def running_lammps():
    import subprocess
    import glob
    import os
    
    # 1. SETUP
    lammps_cmd = "lmp" 
    # Since 'inputfile.txt' is INSIDE the folder, we just use the filename
    input_script = "inputfile.txt" 
    
    # 2. FIND FOLDERS
    # We look for directories named "disp-*"
    # This finds ['disp-001', 'disp-002', ...]
    folders = sorted(glob.glob("disp-*"))
    
    print(f"Found {len(folders)} folders. Starting simulations...")
    
    # 3. RUN LOOP
    for i, folder in enumerate(folders):
        # We assume the file inside matches the folder name
        # e.g., folder "disp-001" contains "disp-001.lammps"
        # We extract the "001" part from the folder name
        idx = folder.split("-")[-1] 
        data_file = f"disp-001.lammps" ##for this run only
        
        # Check if the file actually exists inside the folder
        full_path = os.path.join(folder, data_file)
        if not os.path.exists(full_path):
            print(f"⚠️ Skipping {folder}: {data_file} not found.")
            continue

        # Define log file name (saved inside the folder to keep it clean)
        log_name = f"log.lammps"
        
        # Construct the command
        # Note: We do NOT put "disp-001/" in the filename here.
        # Why? Because we are going to set cwd=folder below.
        cmd = [
            lammps_cmd,
            "-in", input_script,
            "-log", log_name,
            "-var", "fname", data_file
        ]
        
        print(f"Running in {folder} ({i+1}/{len(folders)})...")
        
        try:
            # --- CRITICAL CHANGE ---
            # cwd=folder tells Python: "cd into this folder, then run the command"
            # This fixes all path issues with graph-compress.pb and inputfile.txt
            subprocess.run(cmd, cwd=folder, check=True, text=True) 
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Error in {folder}: {e}")
            # Optional: continue to next file even if one fails
            continue 
    
    print("✅ All runs complete.")

# Run it
if __name__ == "__main__":
    running_lammps()