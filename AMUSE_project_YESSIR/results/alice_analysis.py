#%%
import os
import matplotlib.pyplot as plt
import numpy as np
#%%


# Function to extract data from a text file
def extract_data_from_file(file_path):
    data = np.loadtxt(file_path)  # Assumes data in the file is formatted properly
    mass = data  # Assuming the first column is x-values
    time = np.arange(1, len(data) + 1)
    return mass, time



# extract data from runs with different random seeds
main_directory = "./alice/random_seeds"
folders = os.listdir(main_directory)

path_to_txt = []
all_path = []
for folder in folders:
    folder_path = os.path.join(main_directory, folder)
    if os.path.isdir(folder_path):
        files = os.listdir(folder_path)
        for file in files:
            file_path = os.path.join(folder_path, file)
            all_path.append(file_path)

all_path = np.reshape(all_path,(19,5))

#calculatating average and standard deviation for total accreted mass
#with respect to each velocity
total_accretion = []
error = []
for i in range(5):
    mass_diff = []
    for j in range(19):
        mass,time = extract_data_from_file(all_path[j,i])
        dm = np.sum(mass[-1] - mass[0])
        mass_diff.append(dm)
    mean = np.mean(mass_diff)
    std = np.std(mass_diff)
    total_accretion.append(mean)
    error.append(std)


#Plot result with errors 

fig_path = os.path.join(main_directory,"total mass trend.png")
v = np.array([20,30,40,50,60]) 

plt.errorbar(v,total_accretion,yerr=error, fmt='o', \
             markersize=8, capsize=5, label='Total accretion with Error Bars')

plt.xlabel('Cluster velocity [MSun]]')
plt.ylabel('Total mass accreted [MSun]')
plt.title('Total mass accretion against different velocities')
plt.legend() 
plt.grid(True)
plt.savefig(fig_path)
plt.show()
# %%
