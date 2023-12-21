#%%
import os
os.chdir('../src')

from amuse.units import units

import os
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

from plotters import plot_evolution_mass_accretion, plot_relative_mass

#%%


# Function to extract data from a text file
def extract_data_from_file(file_path):
    data = np.loadtxt(file_path)  # Assumes data in the file is formatted properly
    mass = data  # Assuming the first column is x-values
    time = np.arange(1, len(data) + 1)
    return mass, time

#%%


# extract data from runs with different random seeds
main_directory = "../results/alice/random_seeds"
folders = os.listdir(main_directory)

path_to_txt = []
all_path = []
for folder in folders:
    folder_path = os.path.join(main_directory, folder)
    if os.path.isdir(folder_path):
        files = os.listdir(folder_path)
        files = natsorted(files)
        for file in files:
            file_path = os.path.join(folder_path, file)
            all_path.append(file_path)

all_path = np.reshape(all_path,(19,10))
print(all_path)
#calculatating average and standard deviation for total accreted mass
#with respect to each velocity
total_accretion = []
error = []
for i in range(10):
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
#%%
fig_path = os.path.join(main_directory,"total mass trend.png")
v = np.array([20,25,30,35,40,45,50,55,60,65]) 

#%%
from scipy.optimize import curve_fit

# Examine fitting 2 different curves
def inverse_proportionality(x, a, b):
    
    y = a/x + b
    return y

def direct_proportionality(x, a, b):
    
    y = a*x +b
    return y


inverse_coeff, _ = curve_fit(inverse_proportionality, v, total_accretion)

direct_coeff, _ = curve_fit(direct_proportionality, v, total_accretion)

#%%

x_range = np.linspace(min(v)-5, max(v)+5, num = 100)
# m_accr inversely proportional to v
fit1 = inverse_proportionality(x_range, inverse_coeff[0], inverse_coeff[1])
# m_accr directly proportional (negative) to v
fit2 = direct_proportionality(x_range, direct_coeff[0], direct_coeff[1])


fig, ax = plt.subplots(figsize=(8, 6))
ax.set_facecolor('whitesmoke')
ax.grid(alpha=0.5)
# best fit for inverse proportionality
plt.plot(x_range, fit1, linestyle = '--', color = 'red',
         label = r'Best fit ($\dot{m} \propto$v$^{-1}$)')

# best fit for direct proportionality
#plt.plot(x_range, fit1, linestyle = '--', color = 'green',
#         label = r'Best fit ($\dot{m} \propto$ -v)')

# Simulation data
plt.errorbar(v,total_accretion,yerr=error, fmt='o', \
             markersize=8, capsize=5, color = '#0c2577',
             label='Total accretion with Error Bars')


plt.xlabel('Cluster velocity [km/s]', size = 'x-large')
plt.ylabel('Total mass accreted [MSun]', size = 'x-large')
plt.title('Total mass accretion against different velocities', size = 'x-large')
ax.legend(fontsize = 14) 

plt.tight_layout()
plt.savefig(fig_path)
# %%
collision_velocity = 20
file = '../results/alice/10000_stars/20 kms/Sink mass_20.00000000000005.txt'
sinks_mass_evolution = np.loadtxt(file)
sinks_mass_evolution = sinks_mass_evolution

diff = sinks_mass_evolution[-1]- sinks_mass_evolution[0]

plot_evolution_mass_accretion(sinks_mass_evolution = sinks_mass_evolution, 
                                end_time = 2.8 | units.Myr, 
                                time_step = 0.1 | units.Myr, 
                                velocity = 20)

plot_relative_mass(sinks_mass_evolution = sinks_mass_evolution, 
                    velocity = 20
                    )









