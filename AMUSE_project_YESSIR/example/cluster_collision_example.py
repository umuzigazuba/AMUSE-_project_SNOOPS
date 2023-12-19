# %%

import os
os.chdir('../src')

from amuse.units import units

import socket
from time import time
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from molecular_cloud_initialization import make_molecular_cloud, evolve_molecular_cloud
from utils import make_cluster_with_posvel, code_bridge_channel_initaization, cluster_cloud_collision

# %%

cloud_radius = 15 | units.pc
cloud_density = 10 * 2.3 | (units.amu / units.cm**3)

# Given the cloud radius and density, calculate its total mass and number of SPH particles 
cloud_mass = cloud_density * 4 * np.pi * cloud_radius**3 / 3
cloud_N = cloud_mass.value_in(units.MSun) * 15 # Determined so that one particle has a mass ~ 0.06 MSun
print(f"Cloud with {cloud_mass.in_(units.MSun)}.")
print(f"Number of SPH particles {cloud_N}.")

# %%

initial_cloud, converter_cloud  = make_molecular_cloud(N_cloud = cloud_N,
                                                       M_cloud = cloud_mass,
                                                       R_cloud = cloud_radius,
                                                       seed = 32839)

initial_cloud, density_map = evolve_molecular_cloud(initial_cloud, 
                                                    converter_cloud, 
                                                    end_time = 2 | units.Myr, 
                                                    time_step = 0.2 | units.Myr,  
                                                    seed = 6829)

print("Mean density of the molecular cloud", np.mean(initial_cloud.density))

# %%

velocity = 20
end_time = 2.8 | units.Myr
time_step = 0.1 | units.Myr

print(f"Starting with cluster velocity {velocity}.")

# %%

particles_cloud = initial_cloud.copy()

# %%

star_cluster, converter_cluster = make_cluster_with_posvel(position = 30, # in parsecs
                                                            velocity = velocity, 
                                                            W0 = 7, 
                                                            random_seed = 207349, 
                                                            number_of_stars = 1000)

# %%

sinks, gravity_code, stellar_evolution_code, hydro_code, channels, \
gravity_hydro_bridge  = code_bridge_channel_initaization(time_step, star_cluster, converter_cluster, 
                                                            particles_cloud, converter_cloud, seed = 2804)
# %%

directory_path = '../results/example/'
os.makedirs(directory_path, exist_ok = True)

if __name__ == "__main__":
    start_time = time()
    date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
    print("Collision started on {}".format(socket.gethostname()), f"at time {format(date_time)}")
    
    cluster_post_collision = cluster_cloud_collision(end_time, time_step, sinks, particles_cloud, 
                                                     gravity_code, stellar_evolution_code, hydro_code, 
                                                     channels, gravity_hydro_bridge, directory_path, density_map)

    print("Collision finished (running time: {0:.1f} s)".format(time() - start_time))

# %%
