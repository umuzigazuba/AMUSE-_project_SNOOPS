#%%
import os
os.chdir('../src')

from time import time
from datetime import datetime
import socket
import numpy as np
import matplotlib.pyplot as plt

from molecular_cloud_initialization import *
import utils

from amuse.units import units
#%%

tot_cloud_mass = 4/3 *units.constants.pi * (15 | units.pc)**3 * ( 2.3 | units.amu * 10 / (units.cm**3))
tot_cloud_mass=tot_cloud_mass.value_in(units.MSun)
n_particles = tot_cloud_mass*15
print(tot_cloud_mass,n_particles)
#%%

init_cloud, init_cloud_converter  = make_molecular_cloud(N_cloud = n_particles,
                                                         M_cloud = tot_cloud_mass | units.MSun,
                                                         R_cloud = 15 | units.pc,
                                                         seed = 32839)

init_cloud, density_map = evolve_molecular_cloud(init_cloud, 
                                                    init_cloud_converter, 
                                                    t_end = 1 | units.Myr, 
                                                    dt = 0.2 | units.Myr, 
                                                    seed = 6829)

print("Mean density of the moelcular cloud", np.mean(init_cloud.density))

#%%
velocity = np.array([20])

for i in range (len(velocity)):
    v = velocity[i]
    print("Starting with cluster velocity",v)

    star_cluster,converter_cluster = utils.make_cluster_with_vinit(v,25,207349)

    gravhydrobridge,sinks,channel,particles_cloud,gravity_code,\
        hydro_cloud = utils.AMUSE_bridge_initialization(star_cluster,converter_cluster,\
                                                    init_cloud,init_cloud_converter)

    path = '../results/tests_with_different_criteria/ff_ratio_0.5/'
    variable_name = str(v) + " kms"
    directory_path = os.path.join(path, str(variable_name))
    os.makedirs(directory_path, exist_ok=True)

    if __name__ == "__main__":
        start_time = time()
        date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
        print("Collision started on {}".format(socket.gethostname()),f"at time {format(date_time)}")
        
        post_collision_cluster = utils.let_them_collide_and_save(v,directory_path,2.4,0.2,sinks,gravhydrobridge,channel,\
                                                particles_cloud,gravity_code,hydro_cloud)

        print("Collision finished (running time: {0:.1f}s)".format(time() - start_time))

#%%
#retrieve the data and plot trends 
file_path = '../results/tests_with_different_criteria/ff_ratio_0.5/20 kms/Sink mass_19.977764918756957.txt' 
data = np.loadtxt(file_path)
x_values = np.arange(1, len(data) + 1)

plt.plot(x_values, data)
plt.xlabel('Time Period')
plt.ylabel('Accretion')
plt.title('Accretion Trend')
plt.grid(True)
plt.show()


relative_mass = data - data[0,:]

# Remove the column with the most increase
plt.plot(x_values, relative_mass)
plt.xlabel('Time (0.1 Myr)')
plt.ylabel('Accretion')
plt.grid(True)
plt.show()

# %%
