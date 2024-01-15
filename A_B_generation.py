# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 16:00:26 2024

@author: Khalid
"""

import numpy as np


def AB(interface_along_z,device_length,edge_A,edge_B,reservoir_A,reservoir_B):

    #%% inputs
    
    # interface_along_z=5;
    # device_length=100;
    
    # edge_Al=8;
    # edge_c=8;
    # reservoir_c=15;
    # reservoir_Al=47;
    
    device_length= device_length+reservoir_A+ reservoir_B+edge_A+edge_B
    
    #%%
    
    with open('C_left_interface.txt', 'r') as file:
        # Your code to read or manipulate the file goes here
        content_interface_left = file.read()
    
    lines_interface_left = content_interface_left.splitlines()
    
    int_left=lines_interface_left[interface_along_z]
    
    interface_left_length_z=int_left.split()
    interface_conv_cell=float(interface_left_length_z[1])
    number_atoms_left_interface=lines_interface_left[1]
    number_atoms_left_interface=number_atoms_left_interface.split()
    number_atoms_left_interface=float(number_atoms_left_interface[0])
    number_atoms_left_interface=int(number_atoms_left_interface)
    
    atom_matrix_left_interface=lines_interface_left[14:42+number_atoms_left_interface]
    
    atom_matrix_left_interface=np.genfromtxt(atom_matrix_left_interface)
    
    
    #%%
    with open('C_unitcell.txt', 'r') as file:
        # Your code to read or manipulate the file goes here
        content_A_unitcell = file.read()
    
    lines_A_unitcell = content_A_unitcell.splitlines()
    
    A_unit=lines_A_unitcell[interface_along_z]
    
    A_unit_z=A_unit.split()
    A_unit_cell=float(A_unit_z[1])
    number_atoms_A_unit=lines_A_unitcell[1]
    number_atoms_A_unit=number_atoms_A_unit.split()
    number_atoms_A_unit=float(number_atoms_A_unit[0])
    number_atoms_A_unit=int(number_atoms_A_unit)
    
    atom_matrix_A_unit=lines_A_unitcell[14:42+number_atoms_A_unit]
    
    atom_matrix_A_unit=np.genfromtxt(atom_matrix_A_unit)
    
    #%%
    
    number_A_units=int(np.floor(0.5*device_length/A_unit_cell))
    #A=np.zeros((number_Al_units*atom_matrix_Al_unit.shape[0],5))
    A=atom_matrix_A_unit.copy();
    new_A=atom_matrix_A_unit.copy();
    
    
    for i in range(1,number_A_units-1,1):
        new_A[:,0]=new_A[:,0]+atom_matrix_A_unit.shape[0]
        new_A[:,interface_along_z-1]=new_A[:,interface_along_z-1]+A_unit_cell
        A=np.vstack((A,new_A));
        
    #%%
    
    with open('middle_part.txt', 'r') as file:
        # Your code to read or manipulate the file goes here
        content_middle_part = file.read()
    
    lines_middle_part = content_middle_part.splitlines()
    
    middle_part=lines_middle_part[interface_along_z]
    
    middle_part_z=middle_part.split()
    middle_part_cell=float(middle_part_z[1])
    number_atoms_middle_part=lines_middle_part[1]
    number_atoms_middle_part=number_atoms_middle_part.split()
    number_atoms_middle_part=float(number_atoms_middle_part[0])
    number_atoms_middle_part=int(number_atoms_middle_part)
    
    atom_matrix_middle_part=lines_middle_part[14:42+number_atoms_middle_part]
    
    atom_matrix_middle_part=np.genfromtxt(atom_matrix_middle_part)
    
    
    
    
    
    
    
        
    
    #%%
    
    
    with open('Mo_unitcell.txt', 'r') as file:
        # Your code to read or manipulate the file goes here
        content_B_unitcell = file.read()
    
    lines_B_unitcell = content_B_unitcell.splitlines()
    
    B_unit=lines_B_unitcell[interface_along_z]
    
    B_unit_z=B_unit.split()
    B_unit_cell=float(B_unit_z[1])
    number_atoms_B_unit=lines_B_unitcell[1]
    number_atoms_B_unit=number_atoms_B_unit.split()
    number_atoms_B_unit=float(number_atoms_B_unit[0])
    number_atoms_B_unit=int(number_atoms_B_unit)
    
    atom_matrix_B_unit=lines_B_unitcell[14:42+number_atoms_B_unit]
    
    atom_matrix_B_unit=np.genfromtxt(atom_matrix_B_unit)
    
    #%%
    
    
    number_B_units=int(np.floor(0.5*device_length/B_unit_cell))
    #Al=np.zeros((number_Al_units*atom_matrix_Al_unit.shape[0],5))
    B=atom_matrix_B_unit.copy();
    new_B=atom_matrix_B_unit.copy();
    
    
    for i in range(1,number_B_units-1,1):
        new_B[:,0]=new_B[:,0]+atom_matrix_B_unit.shape[0]
        new_B[:,interface_along_z-1]=new_B[:,interface_along_z-1]+B_unit_cell
        B=np.vstack((B,new_B));
        
    #%%
    
    with open('MO_right_interface.txt', 'r') as file:
        # Your code to read or manipulate the file goes here
        content_interface_right = file.read()
    
    lines_interface_right = content_interface_right.splitlines()
    
    interface_right=lines_interface_right[interface_along_z]
    
    interface_right_z=interface_right.split()
    interface_right_cell=float(interface_right_z[1])
    number_atoms_interface_right=lines_interface_right[1]
    number_atoms_interface_right=number_atoms_interface_right.split()
    number_atoms_interface_right=float(number_atoms_interface_right[0])
    number_atoms_interface_right=int(number_atoms_interface_right)
    
    atom_matrix_interface_right=lines_interface_right[14:42+number_atoms_interface_right]
    
    atom_matrix_interface_right=np.genfromtxt(atom_matrix_interface_right)
    
    
    
    #%%
    
    part1=atom_matrix_left_interface.shape[0]
    A[:,0]=A[:,0]+part1
    A[:,interface_along_z-1]=A[:,interface_along_z-1]+interface_conv_cell
    
    interfce_left_A=np.vstack((atom_matrix_left_interface,A));
    
    part2=interfce_left_A.shape[0]
    
    atom_matrix_middle_part[:,0]=atom_matrix_middle_part[:,0]+part2
    atom_matrix_middle_part[:,interface_along_z-1]=atom_matrix_middle_part[:,interface_along_z-1]+A_unit_cell*(number_A_units-1) +interface_conv_cell
    
    interface_left_A_middle_part=np.vstack((interfce_left_A,atom_matrix_middle_part));
    
    part3=interface_left_A_middle_part.shape[0]
    
    B[:,0]=B[:,0]+part3;
    B[:,interface_along_z-1]=B[:,interface_along_z-1]+ interface_conv_cell+A_unit_cell*(number_A_units-1) +middle_part_cell
    
    interface_left_A_middle_part_B=np.vstack((interface_left_A_middle_part,B));
    
    part4=interface_left_A_middle_part_B.shape[0]
    
    atom_matrix_interface_right[:,0]=atom_matrix_interface_right[:,0]+ part4;
    atom_matrix_interface_right[:,interface_along_z-1]=atom_matrix_interface_right[:,interface_along_z-1]+interface_conv_cell+A_unit_cell*(number_A_units-1) +middle_part_cell +B_unit_cell*(number_B_units-1)
    
    
    AB=np.vstack((interface_left_A_middle_part_B,atom_matrix_interface_right));
    
    
    cell_size_z=interface_conv_cell+A_unit_cell*(number_A_units-1) +middle_part_cell +B_unit_cell*(number_B_units-1)+interface_right_cell
    print(cell_size_z)
    
    
    
    
    if interface_along_z==5:
        
        with open('lammps.txt', 'w') as file:
        # Write content to the file
            file.write('# LAMMPS data file written by OVITO Basic 3.7.12\n')
            file.write('%d atoms\n' % AB.shape[0])
            file.write("%d atom types\n"% int(max(AB[:,1])))
            file.write(lines_A_unitcell[3])
            file.write("\n")
            file.write(lines_A_unitcell[4])
            file.write("\n")
            file.write('0.0 %f zlo zhi \n\n'%cell_size_z)
            file.write('Masses\n \n')
            file.write(lines_A_unitcell[9])
            file.write("\n")
            file.write(lines_A_unitcell[10])
            file.write("\n")
            file.write(lines_A_unitcell[11])
            file.write("\n \n")
            file.write('Atoms  # atomic\n\n')
            for row in AB:
            # Join the row elements with tabs and write to the file
                formatted_row = '\t'.join(['%d' % val if val.is_integer() else '%f' % val for val in row])
                file.write(formatted_row + '\n')
    else:
        with open('lammps.txt', 'w') as file:
        # Write content to the file
            file.write('# LAMMPS data file written by OVITO Basic 3.7.12\n')
            file.write('%d atoms\n' % AB.shape[0])
            file.write("%d atom types\n"% int(max(AB[:,1])))
            file.write('0.0 %f xlo xhi \n'%cell_size_z)
            file.write(lines_A_unitcell[4])
            file.write("\n")
            file.write(lines_A_unitcell[5])
            file.write("\n\n")
            file.write('Masses\n \n')
            file.write(lines_A_unitcell[9])
            file.write("\n")
            file.write(lines_A_unitcell[10])
            file.write("\n")
            file.write(lines_A_unitcell[11])
            file.write("\n \n")
            file.write('Atoms  # atomic\n\n')
            for row in AB:
            # Join the row elements with tabs and write to the file
                formatted_row = '\t'.join(['%d' % val if val.is_integer() else '%f' % val for val in row])
                file.write(formatted_row + '\n')
        
    
    
    return AB, cell_size_z


#%%




#%%

def generate_lammps_input_z(dt, d, T, posi_leftedge, posi_rightedge, posi_hot, posi_cold,nvt,nve):
    input_script = f"""units metal

variable 	T equal {T}
variable 	V equal vol
variable 	dt equal {dt}
variable 	d equal {d}	# thermo interval

# convert from LAMMPS metal units to SI
variable 	kB equal 1.3806504e-23	 # [J/K] Boltzmann constant
variable 	eV2J equal 1.602176565e-19
variable	ps2s equal 1.0e-12
variable 	A2m equal 1.0e-10

# set up problem
dimension 	3
boundary 	p p s
atom_style atomic

neighbor	2.0 nsq
read_data lammps.txt

pair_style mlip mlip.ini
pair_coeff * * 

# define cold and hot regions
variable 	l_c equal 5.442
variable 	Neach equal 80
variable 	Nedge equal 8
variable 	Nres equal 8
variable 	NSED equal 8
variable 	Nabdon equal 4

variable     posi_leftedge equal {posi_leftedge}
variable     posi_rightedge equal {posi_rightedge}
variable     posi_hot equal {posi_hot}
variable     posi_cold equal {posi_cold}

region 112 block INF INF INF INF INF ${{posi_leftedge}} units box
group edge_left region 112
region 113 block INF INF INF INF ${{posi_rightedge}} INF units box
group edge_right region 113
region 114 block INF INF INF INF ${{posi_leftedge}}  ${{posi_rightedge}} units box
group mid region 114

region 2 block INF INF INF INF ${{posi_leftedge}}  ${{posi_hot}} units box
group hot region 2
region 3 block  INF INF INF INF ${{posi_cold}} ${{posi_rightedge}} units box
group cold region 3

velocity edge_left set  0.0 0.0 0.0
velocity edge_right set 0.0 0.0 0.0

# initial velocity
velocity 	mid create 300 728462 mom yes rot yes dist gaussian

timestep 	${{dt}}
thermo_style    custom step temp press pe ke lx ly lz
thermo	${{d}}

# ------------- Equilibration and thermalisation ---------------- 
fix  NVT mid nvt temp ${{T}} ${{T}} 0.10 
dump 1112 all atom 250 nvt.tj #it was npt before
run 	{nvt}
unfix   NVT
undump  1112

# --------------- Equilibration in nve ----------------- 
reset_timestep  0

fix 11 hot  langevin 330.0 330.0 0.5 48279  tally yes
fix 12 cold langevin 270.0 270.0 0.5 48279  tally yes

fix 	NVE mid nve

thermo_style    custom step temp press pe ke lx ly lz f_11 f_12
thermo	${{d}}

run 	{nve}

# output temperature data for post-processing
reset_timestep  0

compute      myKE mid ke/atom    # per atom kinetic energy for temperature calculation
variable temp1 atom c_myKE*${{eV2J}}*2/3/${{kB}}    # 1/2*m*v^2 = 3/2*kB*T, unit conversion

compute cc1 mid chunk/atom bin/1d z lower 2.0 units box
fix 104 mid ave/chunk 10 {int(nve/10)} {nve} cc1 v_temp1 file Tgradient.txt

run {nve}
"""

        

    return input_script

def generate_lammps_input_x(dt, d, T, posi_leftedge, posi_rightedge, posi_hot, posi_cold,nvt,nve):
    input_script = f"""units metal

variable 	T equal {T}
variable 	V equal vol
variable 	dt equal {dt}
variable 	d equal {d}	# thermo interval

# convert from LAMMPS metal units to SI
variable 	kB equal 1.3806504e-23	 # [J/K] Boltzmann constant
variable 	eV2J equal 1.602176565e-19
variable	ps2s equal 1.0e-12
variable 	A2m equal 1.0e-10

# set up problem
dimension 	3
boundary 	s p p
atom_style atomic

neighbor	2.0 nsq
read_data lammps.txt

pair_style mlip mlip.ini
pair_coeff * * 

# define cold and hot regions
variable 	l_c equal 5.442
variable 	Neach equal 80
variable 	Nedge equal 8
variable 	Nres equal 8
variable 	NSED equal 8
variable 	Nabdon equal 4

variable     posi_leftedge equal {posi_leftedge}
variable     posi_rightedge equal {posi_rightedge}
variable     posi_hot equal {posi_hot}
variable     posi_cold equal {posi_cold}

region 112 block INF ${{posi_leftedge}} INF INF INF INF units box
group edge_left region 112
region 113 block ${{posi_rightedge}} INF INF INF INF INF units box
group edge_right region 113
region 114 block ${{posi_leftedge}}  ${{posi_rightedge}} INF INF INF INF units box
group mid region 114

region 2 block ${{posi_leftedge}}  ${{posi_hot}} INF INF INF INF units box
group hot region 2
region 3 block  ${{posi_cold}} ${{posi_rightedge}} INF INF INF INF units box
group cold region 3

velocity edge_left set  0.0 0.0 0.0
velocity edge_right set 0.0 0.0 0.0

# initial velocity
velocity 	mid create 300 728462 mom yes rot yes dist gaussian

timestep 	${{dt}}
thermo_style    custom step temp press pe ke lx ly lz
thermo	${{d}}

# ------------- Equilibration and thermalisation ---------------- 
fix  NVT mid nvt temp ${{T}} ${{T}} 0.10 
dump 1112 all atom 250 nvt.tj #it was npt before
run 	{nvt}
unfix   NVT
undump  1112

# --------------- Equilibration in nve ----------------- 
reset_timestep  0

fix 11 hot  langevin 330.0 330.0 0.5 48279  tally yes
fix 12 cold langevin 270.0 270.0 0.5 48279  tally yes

fix 	NVE mid nve

thermo_style    custom step temp press pe ke lx ly lz f_11 f_12
thermo	${{d}}

run 	{nve}

# output temperature data for post-processing
reset_timestep  0

compute      myKE mid ke/atom    # per atom kinetic energy for temperature calculation
variable temp1 atom c_myKE*${{eV2J}}*2/3/${{kB}}    # 1/2*m*v^2 = 3/2*kB*T, unit conversion

compute cc1 mid chunk/atom bin/1d x lower 2.0 units box
fix 104 mid ave/chunk 10 {int(nve/10)} {nve} cc1 v_temp1 file Tgradient.txt

run {nve}
"""

        

    return input_script



# Example usage


interface_along_z=5;
device_length=500;

edge_A=8;
edge_B=8;
posi_hot=21.0
reservoir_A=posi_hot-edge_A;
reservoir_B=281.0-244.0;


dt = 0.001
d = 2000
[ALC,CELL_SIZE]=AB(interface_along_z, device_length, edge_A, edge_B, reservoir_A, reservoir_B)

posi_rightedge=CELL_SIZE-edge_B;
posi_cold=posi_rightedge-reservoir_B;
T=300;
posi_leftedge=8.0;

nvt=1500000;
nve=1000000;
if interface_along_z==5:
    
    generated_script = generate_lammps_input_z(dt, d, T, posi_leftedge, posi_rightedge, posi_hot, posi_cold,nvt,nve);
else:
    generated_script = generate_lammps_input_x(dt, d, T, posi_leftedge, posi_rightedge, posi_hot, posi_cold,nvt,nve);
with open('inputfile.txt', 'w') as file:
    file.write(generated_script)

    





































