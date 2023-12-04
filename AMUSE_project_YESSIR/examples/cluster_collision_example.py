#%%
import os
os.chdir('../src')

from time import time
from datetime import datetime
import socket
import numpy as np
import matplotlib.pyplot as plt

from molecular_cloud_initialization import *
from plotters import *
from cluster_cloud_initialization import *
from utils import *

from amuse.units import units
#%%


tot_cloud_mass = 4/3 *units.constants.pi * (15 | units.pc)**3 * ( 2.3 | units.amu * 20 / (units.cm**3))
print(tot_cloud_mass.value_in(units.MSun))

#Assmuing a pure molecular hydrogen cloud, with typical density around 85 molecules per cm^-3, calculate the approximate cloud mass based on
#cloud size. 50 pc is selected for a small GC of only 100 stars. 

# %%
# initialise and evolve the MC particle set
init_cloud, init_cloud_converter  = make_molecular_cloud(N_cloud = 100_000,
                                                         M_cloud = 16_000 | units.MSun,
                                                         R_cloud = 15 | units.pc,
                                                         seed = 32839)

init_cloud, density_map = evolve_molecular_cloud(init_cloud, 
                                                    init_cloud_converter, 
                                                    t_end = 2 | units.Myr, 
                                                    dt = 0.2 | units.Myr, 
                                                    seed = 6829)

print("Mean density of the moelcular cloud", np.mean(init_cloud.density))


#%%

velocity = np.array([45])

for i in range (len(velocity)):
    v = velocity[i]
    print("Starting with cluster velocity",v)

    star_cluster,converter_cluster = make_cluster_with_vinit(v,30,207349)

    gravhydrobridge,sinks,channel,particles_cloud,gravity_code,\
        hydro_cloud = AMUSE_bridge_initialization(star_cluster,converter_cluster,\
                                                  init_cloud,init_cloud_converter)
    
    path = '../results/'
    variable_name = v
    directory_path = os.path.join(path, str(variable_name))
    os.makedirs(directory_path, exist_ok=True)

    if __name__ == "__main__":
        start_time = time()
        date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
        print("Collision started on {}".format(socket.gethostname()),f"at time {format(date_time)}")
        
        let_them_collide(directory_path,1.8,0.2,sinks,gravhydrobridge,channel,\
                         particles_cloud,gravity_code,hydro_cloud,density_plot_flag=1)

        print("Collision finished (running time: {0:.1f}s)".format(time() - start_time))
#%%
#retrieve the data and plot trends after removing the most massive star 
file_path = '../results/45/Sink mass_44.97776491875698.txt'  # Replace with your file path

# Read the data from the text file
data = np.loadtxt(file_path)

# Create x-axis values (assuming equal intervals, starting from 1)
x_values = np.arange(1, len(data) + 1)

# Plot the growth trend
plt.plot(x_values, data)
plt.xlabel('Time Period')
plt.ylabel('Accretion')
plt.title('Accretion Trend')
plt.grid(True)
plt.show()

column_increases = np.diff(data, axis=0)
total_increases = np.sum(column_increases, axis=0)

# Find the index of the column with the most increase
column_index_with_max_increase = np.argmax(total_increases)

# Remove the column with the most increase
data_without_max_increase = np.delete(data, column_index_with_max_increase, axis=1)
plt.plot(x_values, data_without_max_increase)
plt.xlabel('Time (0.1 Myr)')
plt.ylabel('Accretion')
plt.grid(True)
plt.show()

# %%

plt.hist(total_increases)
plt.show()
# %%
print(total_increases)
# %%
