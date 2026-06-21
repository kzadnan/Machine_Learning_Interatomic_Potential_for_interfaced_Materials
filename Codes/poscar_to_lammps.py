# -*- coding: utf-8 -*-
"""
Created on Sun Dec  7 13:03:42 2025

@author: Owner
"""

def organize_folders():
    
    import glob
    # 1. Find all files that look like "POSCAR-*"
    # This creates a list like ['POSCAR-0', 'POSCAR-1', ...]
    file_list = glob.glob("POSCAR-*")
    count = len(file_list)
   
    import os
    import shutil
    
    # 1. Define the files that need to be copied into EVERY folder
    # Make sure these files exist in your current directory!
    common_files = ["graph-compress.pb", "inputfile.txt", "job"]
    
    # 2. Define how many folders you need
    total_folders = count  # Change this to your actual number
    
    for i in range(1,total_folders+1):
        # Format the index as 000, 001, etc.
        idx = f"{i:03d}"
        
        # Define folder and file names
        folder_name = f"disp-{idx}"
        disp_file = f"disp-{idx}.lammps"
        
        # --- Step A: Create the folder ---
        # exist_ok=True prevents errors if the folder already exists
        os.makedirs(folder_name, exist_ok=True)
        
        # --- Step B: Move the specific disp file into the folder ---
        # We check if the file exists first to avoid crashing
        if os.path.exists(disp_file):
            # uses shutil.move to put the file inside the new folder
            shutil.move(disp_file, os.path.join(folder_name, disp_file))
            print(f"Moved {disp_file} -> {folder_name}/")
        else:
            print(f"⚠️ Warning: {disp_file} not found, skipping move.")
    
        # --- Step C: Copy the common files (potential, input, job) ---
        for file in common_files:
            if os.path.exists(file):
                shutil.copy(file, folder_name)
            else:
                print(f"❌ Error: Source file '{file}' does not exist!")
    
    print("✅ Organization complete.")