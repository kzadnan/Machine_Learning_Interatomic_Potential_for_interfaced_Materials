# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 17:45:57 2023

@author: Khalid
"""

import numpy as np
import matplotlib.pyplot as plt

with open('Tgradient.txt', 'r') as file:
    # Read the content of the file
    for _ in range(4):
        next(file)      
    Temperature_grad=np.loadtxt(file)
    
Temperature_grad_refined=Temperature_grad[(Temperature_grad[:,3]>0)]


pos_along_sim_domain=Temperature_grad_refined[:,1];
temp_along_sim_domain=Temperature_grad_refined[:,3];

adjusting=3

#%% Region defined
lattice_const=1;
with open('inputfile.txt', 'r') as file:
    # Your code to read or manipulate the file goes here
    content = file.read()
    lines_inputfile = content.splitlines()
    timestep=lines_inputfile[4];
    timestep=timestep.split()
    timestep=float(timestep[3])
    posihot_A=lines_inputfile[37];
    posihot_A=posihot_A.split()
    posihot_A=float(posihot_A[3])
    posicold_B=lines_inputfile[38];
    posicold_B=posicold_B.split()
    posicold_B=float(posicold_B[3])
    
    posiright_B=lines_inputfile[36];
    posiright_B=posiright_B.split()
    posiright_B=float(posiright_B[3])
    
    interface_direction=lines_inputfile[16];
    interface_direction=interface_direction.split()
    if interface_direction[1]=='s':
        interface_dir=1
    elif interface_direction[2]=='s':
        interface_dir=2
    else:
        interface_dir=3
    
edge_A=8
reservoir_A=posihot_A-edge_A
# timestep=0.0005 
edge_B=8;
reservoir_B=posiright_B-posicold_B;

with open('lammps.txt', 'r') as file:
    content = file.read()
    lines_lammps = content.splitlines()
    a=lines_lammps[3]
    a=a.split()
    a=float(a[1])
    
    b=lines_lammps[4]
    b=b.split()
    b=float(b[1])
    
    c=lines_lammps[5]
    c=c.split()
    c=float(c[1])
if interface_dir==1:
    area= b*c
elif interface_dir==2:
    area= a*c
else:
    area=b*a
    


edge_length_A=edge_A*lattice_const;
reservoir_length_A=reservoir_A*lattice_const;
edge_length_B=edge_B*lattice_const;
reservoir_length_B=reservoir_B*lattice_const

x1=edge_length_A+reservoir_length_A;
x3=Temperature_grad[-1,1]-edge_length_B-reservoir_length_B
cell_1=(x1+x3)/2 -adjusting


pos_cell_1=cell_1*lattice_const;


region_1=Temperature_grad_refined[(Temperature_grad_refined[:,1]<=edge_length_A+reservoir_length_A)];
coefficients_region_1=np.polyfit(region_1[:,1],region_1[:,3],1)
[slope_region1,intercept_region1]=coefficients_region_1

region_2=Temperature_grad_refined[(Temperature_grad_refined[:,1]>=edge_length_A+reservoir_length_A) &((Temperature_grad_refined[:,1])<=pos_cell_1-10)];
coefficients_region_2=np.polyfit(region_2[:,1],region_2[:,3],1)
[slope_region2,intercept_region2]=coefficients_region_2


region_3=Temperature_grad_refined[(Temperature_grad_refined[:,1]>=pos_cell_1) &((Temperature_grad_refined[:,1])<=Temperature_grad[-1,1]-edge_length_B-reservoir_length_B)];
coefficients_region_3=np.polyfit(region_3[:,1],region_3[:,3],1)
[slope_region3,intercept_region3]=coefficients_region_3

region_4=Temperature_grad_refined[(Temperature_grad_refined[:,1]>=Temperature_grad[-1,1]-edge_length_B-reservoir_length_B)];
coefficients_region_4=np.polyfit(region_4[:,1],region_4[:,3],1)
[slope_region4,intercept_region4]=coefficients_region_4


plt.figure(1)
plt.xlabel('Position along Sim Domain')
plt.ylabel('Temperature along Sim Domain')
plt.axvline(x=edge_length_A+reservoir_length_A, color='red', linestyle='--',)
plt.axvline(x=pos_cell_1, color='red', linestyle='--',)
plt.axvline(x=Temperature_grad[-1,1]-edge_length_B-reservoir_length_B, color='red', linestyle='--',)
plt.plot(pos_along_sim_domain,temp_along_sim_domain,color='black')
plt.plot(region_1[:,1],slope_region1*region_1[:,1]+intercept_region1,color='pink')
plt.plot(edge_length_A+reservoir_length_A,slope_region1*(edge_length_A+reservoir_length_A)+intercept_region1,marker='s',color='black')

plt.plot(region_2[:,1],slope_region2*region_2[:,1]+intercept_region2,color='magenta')
plt.plot(edge_length_A+reservoir_length_A,slope_region2*(edge_length_A+reservoir_length_A)+intercept_region2,marker='s',color='black')
plt.plot(pos_cell_1,slope_region2*(pos_cell_1)+intercept_region2,marker='s',color='black')

plt.plot(region_3[:,1],slope_region3*region_3[:,1]+intercept_region3,color='blue')
plt.plot(pos_cell_1,slope_region3*(pos_cell_1)+intercept_region3,marker='s',color='black')

plt.plot(region_4[:,1],slope_region4*region_4[:,1]+intercept_region4,color='cyan')
plt.plot(Temperature_grad[-1,1]-edge_length_B-reservoir_length_B,slope_region4*(Temperature_grad[-1,1]-edge_length_B-reservoir_length_B)+intercept_region4,marker='s',color='black')

plt.show()

#%% heat flux calculation


with open('log.txt', 'r') as file:
    content = file.read()
    lines=content.splitlines();
    number_runs=0;
    run_ended=0
    for line in range(int(len(lines))):
        if lines[line]=='Step Temp Press PotEng KinEng Lx Ly Lz f_11 f_12 ' or lines[line]=='Step Temp Press PotEng KinEng Lx Ly Lz ':
            number_runs=number_runs+1
            if number_runs==3:
                skip_lines=line
                number_runs=0
                continue
        if lines[line]=='MPI task timing breakdown:':
            run_ended=run_ended+1
            if run_ended==3:
                last_line=line-6
                
with open('log.txt', 'r') as file:
    for _ in range(skip_lines+1):
        next(file)      
    heatflux=np.loadtxt(file,max_rows=last_line-skip_lines,usecols=range(10))

timesteps=heatflux[:,0]
f11=heatflux[:,8]
f12=heatflux[:,9]

coefficients_heatflux_f11=np.polyfit(timesteps,f11,1)
[slope_f11,intercept_f11]=coefficients_heatflux_f11;

coefficients_heatflux_f12=np.polyfit(timesteps,f12,1)
[slope_f12,intercept_f12]=coefficients_heatflux_f12;

f11_linfitted=slope_f11*timesteps +intercept_f11;
f12_linfitted=slope_f12*timesteps +intercept_f12;


plt.figure(2)
plt.plot(timesteps, f11, label='f11', color='blue')
plt.plot(timesteps, f12, label='f12', color='red')
plt.plot(timesteps, f11_linfitted, linestyle='--', color='blue', label='f11 Linear Fit')
plt.plot(timesteps, f12_linfitted, linestyle='--', color='red', label='f12 Linear Fit')
plt.xlabel('Timesteps')
plt.ylabel('Values')
plt.title('Multiple Lines on the Same Plot')
plt.legend()
plt.show()

ev=1/2*( abs(slope_f11) +abs(slope_f12));
Q=ev* 1/timestep *1e12*1.6E-19;





#%% Ge apparent thermal conductivity 


x1=edge_length_A+reservoir_length_A;
x2=pos_cell_1

x3=Temperature_grad[-1,1]-edge_length_B-reservoir_length_B

T1= x1 * slope_region1 +intercept_region1;
T2=x1 * slope_region2 +intercept_region2;

T3=x2 * slope_region2 +intercept_region2;
T4=x2 * slope_region3 +intercept_region3;

T5=x3 * slope_region3 +intercept_region3;
T6=x3 * slope_region4 +intercept_region4;




L_A=abs(x2-x1);
delta_A_app=T1- T3;
K_app_A= Q*1e20*L_A*1e-10/(area*delta_A_app)


L_B=abs(x3-x2);
delta_B_app=T4- T5; #same as conventional case
K_app_B= Q*1e20*L_B*1e-10/(area*delta_B_app)

delta_temp_jump_first_interface=T3-T4;
TBA_first_interface= Q*1e20/(area*delta_temp_jump_first_interface*1e6)

device_length=L_A+L_B

  
    


print("Apparent Thermal conductivity of A is %f W/m-K" %K_app_A);
print("Apparent or conventional Thermal conductivity of B is %f W/m-K" %K_app_B);
print("Thermal boundary conductance of the first interface is %f MW/m^2-k"%TBA_first_interface)
print("Device length of the interface is %f"%device_length)
























