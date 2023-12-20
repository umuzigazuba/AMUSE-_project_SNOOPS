# %%

import os
import matplotlib.pyplot as plt
import numpy as np
from plotters import plot_evolution_mass_accretion, plot_relative_mass
from amuse.units import units
from natsort import natsorted

# %%

# Function to extract data from a text file
def extract_data_from_file(file_path):
    data = np.loadtxt(file_path)  # Assumes data in the file is formatted properly
    mass = data  # Assuming the first column is x-values
    time = np.arange(1, len(data) + 1)
    return mass, time

# %%

# Extract data from runs with different random seeds
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

all_path = np.reshape(all_path, (19, 10))
# print(all_path) # Check that all runs are used to create the plot

# Calculatating average and standard deviation for total accreted mass
# With respect to each velocity
total_accretion = []
error = []
for i in range(10):
    mass_diff = []
    for j in range(19):
        mass, time = extract_data_from_file(all_path[j, i])
        dm = np.sum(mass[-1] - mass[0])
        mass_diff.append(dm)
    mean = np.mean(mass_diff)
    std = np.std(mass_diff)
    total_accretion.append(mean)
    error.append(std)

# Plot result with errors 
fig_path = os.path.join(main_directory,"total mass trend.png")
cluster_velocitiies = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60, 65]) 

fig, ax = plt.subplots(figsize = (7, 5))
ax.set_facecolor('whitesmoke')

plt.errorbar(cluster_velocitiies, total_accretion, yerr = error, fmt = 'o', color = '#0c2577', \
             markersize = 8, capsize = 5, label = 'Total accretion with error bars')

coefficients = np.polyfit(cluster_velocitiies, total_accretion, 1)
poly_function = np.poly1d(coefficients)
x_fit = np.linspace(min(cluster_velocitiies), max(cluster_velocitiies), 100)
plt.plot(x_fit, poly_function(x_fit), color='red',linestyle = '--', label='Fitted Line')

plt.xticks(cluster_velocitiies)
plt.xlabel('Cluster velocity [km/s]')
plt.ylabel('Total mass accreted [MSun]')
plt.title('Total mass accretion against different cluster velocities')
plt.legend() 
plt.grid(True)
plt.savefig(fig_path)
plt.show()
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
                    velocity = 20)

# %%