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
n_particles = tot_cloud_mass*25
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
write_to_file(init_cloud,"molecular_cloud.hdf5")
#%%
import importlib
importlib.reload(utils)

#%%
velocity = np.array([20])

for i in range (len(velocity)):
    v = velocity[i]
    print("Starting with cluster velocity",v)

    star_cluster,converter_cluster = utils.make_cluster_with_vinit(v,30,207349)

    gravhydrobridge,sinks,channel,particles_cloud,gravity_code,\
        hydro_cloud = utils.AMUSE_bridge_initialization(star_cluster,converter_cluster,\
                                                    init_cloud,init_cloud_converter)

    path = '../results/tests_with_tff/0/'
    variable_name = str(v) + " kms"
    directory_path = os.path.join(path, str(variable_name))
    os.makedirs(directory_path, exist_ok=True)

    if __name__ == "__main__":
        start_time = time()
        date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
        print("Collision started on {}".format(socket.gethostname()),f"at time {format(date_time)}")
        
        post_collision_cluster = utils.let_them_collide_and_save(directory_path,2.5,0.1,sinks,gravhydrobridge,channel,\
                                                particles_cloud,gravity_code,hydro_cloud)

        print("Collision finished (running time: {0:.1f}s)".format(time() - start_time))
# %%

