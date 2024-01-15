# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 22:37:29 2024

@author: Khalid
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 16:59:31 2024

@author: Khalid
"""

import numpy as np
import matplotlib.pyplot as plt
with open('total.cfg', 'r') as file:
    # Your code to read or manipulate the file goes here
    content = file.read()

lines = content.splitlines()
size=int(lines[2])

if lines[size+13]=='END_CFG':
    cfg_size=size+15;
else:
    cfg_size=size+16
    
# Print the first 5 lines as an example

size_lines=len(lines)

train_data=0.8
test_data=1-train_data

n_train=int(train_data*size_lines/cfg_size);
n_test=int(test_data*size_lines/cfg_size);

train=[];
test=[];
energy=[];
en_step=size+9;
#%%
for k in range(en_step,size_lines-5,cfg_size):
    energy.append(float(lines[k]))
    
cfgs=np.linspace(1,len(energy),len(energy));
plt.figure(1)
plt.title(" Total cfg energy")
plt.plot(cfgs,(energy),color='r',linestyle='-')


lower=min(energy)
upper=max(energy)

i=0;
for j in range(int(n_train)):
    if float(lines[j*cfg_size+en_step])>=upper or float(lines[j*cfg_size+en_step])<=lower:
        i=i+cfg_size-1;
    else:
        for i in range(int(j*cfg_size),int((j+1)*cfg_size),1):
            train.append(lines[i])        
    

    
energy_train=[];
for kk in range(en_step,len(train),cfg_size):
    energy_train.append(float(train[kk]))    
 
train_cfgs=np.linspace(1,len(energy_train),len(energy_train));

plt.figure(2)
plt.plot(train_cfgs,energy_train,color='b',linestyle='--',label="train") 
plt.title("Energy in train and test database")

for j in range(int(cfg_size*n_train),int(size_lines),1):
    test.append(lines[j])
    
energy_test=[];
for kk in range(en_step,len(test),cfg_size):
    energy_test.append(float(test[kk])) 
    
test_cfgs=np.linspace(1,len(energy_test),len(energy_test)); 
test_cfgs=test_cfgs+n_train   

plt.plot(test_cfgs,energy_test,color='g',linestyle='--',label="test")   
plt.legend()



 
#%%
# Specify the file path
file_path = "train.cfg"

# Open the file and write the list elements
with open(file_path, 'w') as file:
    for item in train:
        file.write(f"{item}\n")

print(f"List elements written to {file_path}")


# Specify the file path
file_path = "test.cfg"

# Open the file and write the list elements
with open(file_path, 'w') as file:
    for item in test:
        file.write(f"{item}\n")

print(f"List elements written to {file_path}")
    
    