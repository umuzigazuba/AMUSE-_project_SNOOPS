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
tot_cloud_mass=tot_cloud_mass.value_in(units.MSun)
n_particles = tot_cloud_mass*20

init_cloud, init_cloud_converter  = make_molecular_cloud(N_cloud = n_particles,
                                                         M_cloud = tot_cloud_mass | units.MSun,
                                                         R_cloud = 15 | units.pc,
                                                         seed = 32839)

init_cloud, density_map = evolve_molecular_cloud(init_cloud, 
                                                    init_cloud_converter, 
                                                    t_end = 2 | units.Myr, 
                                                    dt = 0.2 | units.Myr, 
                                                    seed = 6829)

print("Mean density of the moelcular cloud", np.mean(init_cloud.density))
#%%

#%%

velocity = np.array([50])

# for i in range (len(velocity)):
#   v = velocity[i]
v=50
print("Starting with cluster velocity",v)

star_cluster,converter_cluster = make_cluster_with_vinit(v,30,207349)

gravhydrobridge,sinks,channel,particles_cloud,gravity_code,\
    hydro_cloud = AMUSE_bridge_initialization(star_cluster,converter_cluster,\
                                                init_cloud,init_cloud_converter)

path = '../results/dt=0.2'
variable_name = v
directory_path = os.path.join(path, str(variable_name))
os.makedirs(directory_path, exist_ok=True)

if __name__ == "__main__":
    start_time = time()
    date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
    print("Collision started on {}".format(socket.gethostname()),f"at time {format(date_time)}")
    
    post_collision_cluster,cloud_density_cubes,\
        star_position,xgrid, ygrid, zgrid = let_them_collide_and_save\
                                            (directory_path,1.8,0.2,sinks,gravhydrobridge,channel,\
                                            particles_cloud,gravity_code,hydro_cloud)

    print("Collision finished (running time: {0:.1f}s)".format(time() - start_time))
#%%
#retrieve the data and plot trends 
# file_path = '../results/45/Sink mass_44.97776491875698.txt' 
# data = np.loadtxt(file_path)
# x_values = np.arange(1, len(data) + 1)

# plt.plot(x_values, data)
# plt.xlabel('Time Period')
# plt.ylabel('Accretion')
# plt.title('Accretion Trend')
# plt.grid(True)
# plt.show()

# column_increases = np.diff(data, axis=0)
# total_increases = np.sum(column_increases, axis=0)

# # Find the index of the column with the most increase
# column_index_with_max_increase = np.argmax(total_increases)

# # Remove the column with the most increase
# data_without_max_increase = np.delete(data, column_index_with_max_increase, axis=1)
# plt.plot(x_values, data_without_max_increase)
# plt.xlabel('Time (0.1 Myr)')
# plt.ylabel('Accretion')
# plt.grid(True)
# plt.show()

# diff = (data[-1,:] - data[0,:])
# plt.hist(diff)
# plt.show()

# diff = (data_without_max_increase[-1,:] - data_without_max_increase[0,:])
# plt.hist(diff)
# plt.show()
# %%
