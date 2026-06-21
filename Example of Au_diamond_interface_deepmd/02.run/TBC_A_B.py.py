import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================
# Update these parameters for new LAMMPS runs instead of relying on inputfile.txt
CONFIG = {
    'timestep': 0.0005,         # ps (metal units)
    'lattice_const': 1.0,       # Scaling factor if needed
    'cross_section_area': 7.6233662 * 15.246732, # Lx * Ly from log.txt
    
    # Geometry (in Angstroms)
    'edge_A': 5.0,
    'edge_B': 5.0,
    'posihot_A': 115.17,        # Example coordinate
    'posicold_B': 40.44,        # Example coordinate
    'posiright_B': 125.17,      # End of simulation box
    
    # Analysis adjustments
    'nve_stable_steps': 0,      # Steps to skip at the beginning of production
    'adjusting': 5.0,          # Interface cell adjustment
}

# Derived geometric lengths
reservoir_A = CONFIG['posicold_B'] - CONFIG['edge_A']
reservoir_B = CONFIG['posiright_B'] - CONFIG['posihot_A']

EDGE_A_LEN = CONFIG['edge_A'] * CONFIG['lattice_const']
RES_A_LEN = reservoir_A * CONFIG['lattice_const']
EDGE_B_LEN = CONFIG['edge_B'] * CONFIG['lattice_const']
RES_B_LEN = reservoir_B * CONFIG['lattice_const']


# =============================================================================
# 2. DATA EXTRACTION FUNCTIONS
# =============================================================================
def read_temperature_gradient(filename):
    """Reads chunk-averaged temperature data, skipping headers and zero-counts."""
    # numpy loadtxt can automatically skip the 4 header lines
    temp_data = np.loadtxt(filename, skiprows=4)
    # Filter out bins where atom count is 0
    temp_data_refined = temp_data[temp_data[:, 3] > 0]
    
    pos = temp_data_refined[:, 1]
    temp = temp_data_refined[:, 3]
    return temp_data_refined, pos, temp

def read_lammps_log(filename):
    """Robustly extracts the last thermodynamic block (production run) from log.txt."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    start_idx = -1
    end_idx = -1
    header = []
    
    # Iterate backwards to find the final production run
    for i in range(len(lines)-1, -1, -1):
        if "Loop time" in lines[i]:
            end_idx = i
        if "Step Temp Press" in lines[i]:
            start_idx = i
            header = lines[i].split()
            break
            
    if start_idx == -1 or end_idx == -1:
        raise ValueError("Could not find complete thermodynamic data block in log file.")

    # Parse data into float array
    data = []
    for line in lines[start_idx + 1 : end_idx]:
        data.append([float(x) for x in line.split()])
        
    return np.array(data), header


# =============================================================================
# 3. PHYSICS CALCULATIONS
# =============================================================================
def calculate_heat_flux(thermo_data, header, timestep):
    """Calculates the total heat flux Q (Watts) from cumulative energy tallies."""
    step_col = header.index('Step')
    f11_col = header.index('f_11')
    f12_col = header.index('f_12')
    
    timesteps = thermo_data[:, step_col]
    f11 = thermo_data[:, f11_col]
    f12 = thermo_data[:, f12_col]
    
    # Linear fit to find dE/dt
    slope_f11, int_f11 = np.polyfit(timesteps, f11, 1)
    slope_f12, int_f12 = np.polyfit(timesteps, f12, 1)
    
    # Calculate energy transfer rate (eV/step -> Watts)
    # 1.6e-19 J/eV, 1e12 ps/s
    ev_per_step = 0.5 * (abs(slope_f11) + abs(slope_f12))
    #ev_per_step=abs(slope_f11)
    #ev_per_step=abs(slope_f12)

    
    Q_watts = ev_per_step * (1 / timestep) * 1e12 * 1.6e-19
    
    fits = {
        't': timesteps,
        'f11': f11, 'f11_fit': slope_f11 * timesteps + int_f11,
        'f12': f12, 'f12_fit': slope_f12 * timesteps + int_f12
    }
    
    return Q_watts, fits


# =============================================================================
# 4. MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    
    # --- 1. Load Data ---
    temp_data, pos, temp = read_temperature_gradient('Tgradient.txt')
    thermo_data, thermo_header = read_lammps_log('log.txt')
    
    # --- 2. Define Boundaries ---
    x_left_res_end = EDGE_A_LEN + RES_A_LEN
    x_right_res_start = pos[-1] - EDGE_B_LEN - RES_B_LEN
    pos_interface = ((x_left_res_end + x_right_res_start) / 2) - CONFIG['adjusting']
    
    # --- 3. Piecewise Linear Fitting ---
    masks = [
        (pos <= x_left_res_end),
        (pos >= x_left_res_end) & (pos <= pos_interface - 2),
        (pos >= pos_interface) & (pos <= x_right_res_start),
        (pos >= x_right_res_start)
    ]
    
    regions = []
    for m in masks:
        region_data = temp_data[m]
        slope, intercept = np.polyfit(region_data[:, 1], region_data[:, 3], 1)
        regions.append({'data': region_data, 'slope': slope, 'intercept': intercept})
        
    # --- 4. Calculate Heat Flux ---
    # Strip equilibration steps if defined
    stable_data = thermo_data[CONFIG['nve_stable_steps']:, :]
    Q, flux_fits = calculate_heat_flux(stable_data, thermo_header, CONFIG['timestep'])
    
    # --- 5. Calculate Thermal Properties ---
    # Extrapolate temperatures to boundaries
    T1 = x_left_res_end * regions[0]['slope'] + regions[0]['intercept']
    T2 = x_left_res_end * regions[1]['slope'] + regions[1]['intercept']
    
    T3 = pos_interface * regions[1]['slope'] + regions[1]['intercept']
    T4 = pos_interface * regions[2]['slope'] + regions[2]['intercept']
    
    T5 = x_right_res_start * regions[2]['slope'] + regions[2]['intercept']
    T6 = x_right_res_start * regions[3]['slope'] + regions[3]['intercept']
    
    # Area conversion: Angstrom^2 to m^2
    area_m2 = CONFIG['cross_section_area'] * 1e-20
    
    # Material A
    L_A_m = abs(pos_interface - x_left_res_end) * 1e-10
    delta_A = T1 - T3
    K_app_A = (Q * L_A_m) / (area_m2 * delta_A)
    
    # Material B
    L_B_m = abs(x_right_res_start - pos_interface) * 1e-10
    delta_B = T4 - T5
    K_app_B = (Q * L_B_m) / (area_m2 * delta_B)
    
    # Interface
    delta_T_interface = T3 - T4
    TBC_W = Q / (area_m2 * delta_T_interface) # W/m^2-K
    TBC_MW = TBC_W / 1e6                      # MW/m^2-K
    
    device_length = abs(pos_interface - x_left_res_end) + abs(x_right_res_start - pos_interface)

    # --- 6. Print Results ---
    print(f"Heat Flux (Q): {Q:.2e} Watts")
    print(f"Apparent Thermal Conductivity (Material A): {K_app_A:.3f} W/m-K")
    print(f"Apparent Thermal Conductivity (Material B): {K_app_B:.3f} W/m-K")
    print(f"Thermal Boundary Conductance (Interface):   {TBC_MW:.3f} MW/m^2-K")
    print(f"Total Device Length Analyzed:               {device_length:.3f} Angstroms")

    # --- 7. Plotting ---
    # Temperature Profile
    plt.figure(figsize=(10, 6))
    plt.plot(pos, temp, color='black', alpha=0.5, label='Raw Data')
    colors = ['pink', 'magenta', 'blue', 'cyan']
    
    for i, r in enumerate(regions):
        plt.plot(r['data'][:, 1], r['slope'] * r['data'][:, 1] + r['intercept'], color=colors[i], linewidth=2)
        
    plt.axvline(x=x_left_res_end, color='red', linestyle='--', label='Baths')
    plt.axvline(x=pos_interface, color='green', linestyle='--', label='Interface')
    plt.axvline(x=x_right_res_start, color='red', linestyle='--')
    
    plt.xlabel('Position along Z (Angstroms)')
    plt.ylabel('Temperature (K)')
    plt.title('Piecewise Linear Fit of Temperature Gradient')
    plt.legend()
    plt.show()

    # Heat Flux Profile
    plt.figure(figsize=(10, 6))
    plt.plot(flux_fits['t'], flux_fits['f11'], label='Hot Bath (f11)', color='red', alpha=0.6)
    plt.plot(flux_fits['t'], flux_fits['f12'], label='Cold Bath (f12)', color='blue', alpha=0.6)
    plt.plot(flux_fits['t'], flux_fits['f11_fit'], 'k--', label='Linear Fit')
    plt.plot(flux_fits['t'], flux_fits['f12_fit'], 'k--')
    
    plt.xlabel('Timestep')
    plt.ylabel('Cumulative Energy (eV)')
    plt.title('Cumulative Energy Transfer vs Time')
    plt.legend()
    plt.show()