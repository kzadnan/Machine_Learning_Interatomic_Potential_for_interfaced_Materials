# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 23:37:20 2024

@author: Khalid
"""
import numpy as np
import matplotlib.pyplot as plt

with open('test.cfg', 'r') as file:
    # Your code to read or manipulate the file goes here
    content_dft = file.read()
    
lines_dft = content_dft.splitlines()
size_dft=int(lines_dft[2])
if lines_dft[size_dft+13]=='END_CFG':
    cfg_size=size_dft+15;
else:
    cfg_size=size_dft+16

# Print the first 5 lines as an example

size_lines_dft=len(lines_dft)
en_step=size_dft+9;
energy_dft=[]
for k in range(en_step,size_lines_dft-5,cfg_size):
    energy_dft.append(float(lines_dft[k]))
    
cfgs_dft=np.linspace(1,len(energy_dft),len(energy_dft));
plt.figure(1)
plt.title("DFT and MTP energy vs cfgs")
plt.xlabel("cfgs")
plt.ylabel("Energy(ev)")
plt.plot(cfgs_dft,(energy_dft),color='r',linestyle='--')

#%%

with open('out.cfg', 'r') as file:
    # Your code to read or manipulate the file goes here
    content_mtp = file.read()
    
lines_mtp = content_mtp.splitlines()
size_mtp=int(lines_mtp[2])
cfg_size_mtp=cfg_size 
# Print the first 5 lines as an example

size_lines_mtp=len(lines_mtp)
en_step=size_mtp+9;
energy_mtp=[]
for k in range(en_step,size_lines_mtp-5,cfg_size_mtp):
    energy_mtp.append(float(lines_mtp[k]))
    
cfgs_mtp=np.linspace(1,len(energy_mtp),len(energy_mtp));
plt.plot(cfgs_mtp,(energy_mtp),color='g',linestyle=':')

#%%
plt.figure(2)
plt.title("DFT and MTP energy")
plt.xlabel("DFT energy (ev)")
plt.ylabel("MTP energy (ev)")
plt.plot(energy_dft,energy_mtp,'*')
plt.plot(energy_dft,energy_dft,linestyle='solid',color='y')
en_x=np.linspace(min(energy_dft),max(energy_dft))
plt.plot(en_x,en_x,color='k',linestyle='-')
plt.xlim(min(energy_dft),max(energy_dft))
plt.ylim(min(energy_dft),max(energy_dft))






#%%
stress_step=en_step+2

sigma_xx_dft=[]
sigma_yy_dft=[]
sigma_zz_dft=[]
sigma_yz_dft=[]
sigma_xz_dft=[]
sigma_xy_dft=[]
stress_dft=[]
stress_mtp=[]
for k in range(stress_step,size_lines_dft-3,cfg_size):
    stress_dft.append(lines_dft[k])

stress_dft=np.genfromtxt(stress_dft)

sigma_xx_dft=stress_dft[:,0]
sigma_yy_dft=stress_dft[:,1]
sigma_zz_dft=stress_dft[:,2]
sigma_yz_dft=stress_dft[:,3]
sigma_xz_dft=stress_dft[:,4]
sigma_xy_dft=stress_dft[:,5]

for k in range(stress_step,size_lines_dft-3,cfg_size):
    stress_mtp.append(lines_mtp[k])

stress_mtp=np.genfromtxt(stress_mtp)

sigma_xx_mtp=stress_mtp[:,0]
sigma_yy_mtp=stress_mtp[:,1]
sigma_zz_mtp=stress_mtp[:,2]
sigma_yz_mtp=stress_mtp[:,3]
sigma_xz_mtp=stress_mtp[:,4]
sigma_xy_mtp=stress_mtp[:,5]

plt.figure(3)
plt.title("DFT and MTP stresses")
plt.xlabel(" DFT stress")
plt.ylabel(" MTP stress")

plt.plot(sigma_xx_dft,sigma_xx_mtp,'*')
plt.plot(sigma_yy_dft,sigma_yy_mtp,'*')
plt.plot(sigma_zz_dft,sigma_zz_mtp,'*')
plt.plot(sigma_yz_dft,sigma_yz_mtp,'*')
plt.plot(sigma_xz_dft,sigma_xz_mtp,'*')
plt.plot(sigma_xy_dft,sigma_xy_mtp,'*')
plt.plot(sigma_xx_dft,sigma_xx_dft,color='k')
stress_x=np.linspace(np.min(stress_dft),np.max(stress_dft))
plt.plot(stress_x,stress_x,color='k',linestyle='-')
plt.xlim(np.min(stress_dft),np.max(stress_dft))
plt.ylim(np.min(stress_dft),np.max(stress_dft))

#%% force calculation


fx_dft=[]
fy_dft=[]
fz_dft=[]
fx_mtp=[]
fy_mtp=[]
fz_mtp=[]

plt.figure(4)
plt.xlabel('DFT_forces')
plt.ylabel('MTP_forces')
plt.title("DFT and MTP forces")
for line in range(int(len(lines_dft))):
    if lines_dft[line]==' AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz':
        
        forces_cfg_dft=lines_dft[line+1:line+1+size_dft]
        forces_cfg_dft=np.genfromtxt(forces_cfg_dft)
        
        fx_dft.append(forces_cfg_dft[:,5])
        fy_dft.append(forces_cfg_dft[:,6])
        fz_dft.append(forces_cfg_dft[:,7])
        
        forces_cfg_mtp=lines_mtp[line+1:line+1+size_dft]
        forces_cfg_mtp=np.genfromtxt(forces_cfg_mtp)
        
        fx_mtp.append(forces_cfg_mtp[:,5])
        fy_mtp.append(forces_cfg_mtp[:,6])
        fz_mtp.append(forces_cfg_mtp[:,7])
        
        
forces=np.hstack((fx_dft,fy_dft,fz_dft))
        
plt.plot(fx_dft,fx_mtp,"*")

plt.plot(fy_dft,fy_mtp,"*")

plt.plot(fz_dft,fz_mtp,"*")
plt.plot(fx_dft,fx_dft,"*")

forces_x=np.linspace(np.min(forces),np.max(forces))
plt.plot(forces_x,forces_x,color='k',linestyle='-')
plt.xlim(np.min(forces_x),np.max(forces_x))
plt.ylim(np.min(forces_x),np.max(forces_x))


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

