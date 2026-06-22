import numpy as np
import matplotlib.pyplot as plt

def calculate_rmse(true, pred):
    return np.sqrt(np.mean((true - pred)**2))

# --- 1. Energy Analysis ---
data_e = np.genfromtxt("results.e.out")
dft_energy = data_e[:, 0]
deepmd_energy = data_e[:, 1]
rmse_e = calculate_rmse(dft_energy, deepmd_energy)

plt.figure(1)
plt.title(f"Energy Correlation (RMSE: {rmse_e:.4f} eV)")
plt.xlabel("DFT energy (eV)")
plt.ylabel("DeePMD energy (eV)")
plt.scatter(dft_energy, deepmd_energy, alpha=0.5)
plt.plot(dft_energy, dft_energy, color='k', linestyle='--')
plt.show()

# --- 2. Forces Analysis ---
data_f = np.genfromtxt("results.f.out")
# Columns 0,1,2 are DFT (x,y,z); Columns 3,4,5 are DeePMD (x,y,z)
dft_forces = data_f[:, 0:3]
deepmd_forces = data_f[:, 3:6]
rmse_f = calculate_rmse(dft_forces, deepmd_forces)

plt.figure(2)
plt.title(f"Force Correlation (RMSE: {rmse_f:.4f} eV/Å)")
plt.xlabel("DFT forces")
plt.ylabel("DeePMD forces")
for i in range(3):
    plt.scatter(data_f[:, i], data_f[:, 3+i], alpha=0.3)
plt.plot(data_f[:, 0], data_f[:, 0], color='k', linestyle='--')
plt.show()

# --- 3. Stress Analysis ---
data_v = np.genfromtxt("results.v.out")
# First 9 columns are DFT, next 9 are DeePMD
dft_stress = data_v[:, 0:9]
deepmd_stress = data_v[:, 9:18]
rmse_v = calculate_rmse(dft_stress, deepmd_stress)

plt.figure(3)
plt.title(f"Stress Correlation (RMSE: {rmse_v:.4f} eV)")
plt.xlabel("DFT stress")
plt.ylabel("DeePMD stress")
for i in range(9):
    plt.scatter(data_v[:, i], data_v[:, 9+i], alpha=0.3)
plt.plot(data_v[:, 0], data_v[:, 0], color='k', linestyle='--')
plt.show()

# --- Summary Output ---
print("-" * 30)
print(f"RMSE Energy: {rmse_e:.6f} eV")
print(f"RMSE Forces: {rmse_f:.6f} eV/Å")
print(f"RMSE Stress: {rmse_v:.6f} eV")
print("-" * 30)