# %%

from amuse.units import units

import os
import socket
from time import time
from datetime import datetime

import numpy as np

from molecular_cloud_initialization import make_molecular_cloud, evolve_molecular_cloud
from utils import make_cluster_with_posvel, code_bridge_channel_initaization, cluster_cloud_collision

# %%

cloud_radius = 15 | units.pc # After evolution the radius will be ~ 20 pc
cloud_density = 10 * 2.3 | (units.amu / units.cm**3)

# Given the cloud radius and density, calculate its total mass and number of
# SPH particles 
cloud_mass = cloud_density * 4 * np.pi * cloud_radius**3 / 3
cloud_N = cloud_mass.value_in(units.MSun) * 15 # Determined so that one
                                            # particle has a mass ~ 0.06 MSun
print(f"Cloud with {cloud_mass.in_(units.MSun)}.")
print(f"Number of SPH particles {cloud_N}.")

# %%
# create the molecular cloud (and its converter)
initial_cloud, converter_cloud  = make_molecular_cloud(N_cloud = cloud_N,
                                                       M_cloud = cloud_mass,
                                                       R_cloud = cloud_radius,
                                                       seed = 32839)
# pre-evolve the cloud for 2 Myr
initial_cloud, density_map = evolve_molecular_cloud(initial_cloud, 
                                                    converter_cloud, 
                                                    end_time = 2 | units.Myr, 
                                                    time_step = 0.2 | units.Myr, 
                                                    seed = 6829)

print("Mean density of the molecular cloud", np.mean(initial_cloud.density))

#%%
# velocities and timesteps for consecutive simulation runs
cluster_velocities = [60, 50, 40, 30, 20] # in km/s
end_times = [1.0, 1.2, 1.6, 2.0, 2.8] # in Myr
time_step = 0.1 | units.Myr

final_clusters = []

for idx in range(len(cluster_velocities)):

    velocity = cluster_velocities[idx]
    end_time = end_times[idx] | units.Myr
  
    print(f"Starting with cluster velocity {velocity}.")

    particles_cloud = initial_cloud.copy()
    # pre-evolve a cluster and assign it the position and velocity for the run
    star_cluster, converter_cluster = make_cluster_with_posvel(position = 30, # in parsecs
                                                               velocity = velocity, 
                                                               random_seed = 207349, 
                                                               number_of_stars = 1000)
    # initialise the codes and the bridge
    sinks, gravity_code, stellar_evolution_code, hydro_code, channels, \
    gravity_hydro_bridge  = code_bridge_channel_initaization(time_step, star_cluster, converter_cluster, 
                                                             particles_cloud, converter_cloud, seed = 2804)

    path = '../results/final_with_stellar_evolution'
    variable_name = str(velocity) + " kms"
    directory_path = os.path.join(path, str(variable_name))
    os.makedirs(directory_path, exist_ok = True)

    if __name__ == "__main__":
        start_time = time()
        date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
        print("Collision started on {}".format(socket.gethostname()), f"at time {format(date_time)}")
        # simulation run
        cluster_post_collision = cluster_cloud_collision(end_time, time_step, sinks, particles_cloud, 
                                                         gravity_code, stellar_evolution_code, hydro_code, 
                                                         channels, gravity_hydro_bridge, directory_path, density_map)

        print("Collision finished (running time: {0:.1f} s)".format(time() - start_time))
        final_clusters.append(cluster_post_collision)

# %%
