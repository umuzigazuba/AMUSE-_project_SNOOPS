# %%

import os
os.chdir('../src')

from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.seba.interface import SeBa

import socket
from time import time
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from molecular_cloud_initialization import make_molecular_cloud, evolve_molecular_cloud
from collision_utils import make_cluster_with_posvel, code_bridge_channel_initaization, cluster_cloud_collision
from analysis_utils import updated_metallicity, evolve_single_star
from plotters import metallicity_histogram, HR_diagramme

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
                                                            W0 = 3, 
                                                            random_seed = 207349, 
                                                            number_of_stars = 1000)

# %%

sinks, gravity_code, stellar_evolution_code, hydro_code, channels, \
gravity_hydro_bridge  = code_bridge_channel_initaization(time_step, star_cluster, converter_cluster, 
                                                            particles_cloud, converter_cloud, seed = 2804)
# %%

directory_path = '../example/results_example/'
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

mass_file = 'Sink_mass_20.0.txt'
sinks_mass_snapshots = np.loadtxt(directory_path + mass_file)

initial_masses = sinks_mass_snapshots[0]
final_masses = sinks_mass_snapshots[-1]
# Compute how much mass is acreeted from the molecular cloud
accreted_mass = np.subtract(final_masses, initial_masses)
# Compute the new metallicities of all stars in the cluster
new_metallicities = updated_metallicity(M_star = initial_masses, M_material = accreted_mass)

# Mask the elements that show change in their metallicity
mask = np.where(new_metallicities > 0.002)
altered_metallicities = new_metallicities[mask[0]]


metallicity_histogram(altered_metallicities, N_bins = 30, collision_velocity = velocity, save_to = directory_path)

# %%

# Make the mask for the non-accreting stars
mask_non_accretion = np.where(new_metallicities == 0.002)
# Get the initial masses of the non accreting stars
non_accreted_masses = initial_masses[mask_non_accretion[0]] | units.MSun
# Make a particle set of the old population
old_population = Particles(mass = non_accreted_masses)

# Stellar evolution of old population
stellar_II = SeBa()
stellar_II.parameters.metallicity = 0.002
stellar_II.particles.add_particles(old_population)
channel_II = stellar_II.particles.new_channel_to(old_population)
channel_II.copy()

# Evovle for the age of the cluster + the time spend in the collision with the molecular cloud
end_time = (10 | units.Gyr) + end_time
model_time = 0.0 | units.Gyr
step = 2 | units.Gyr

while model_time < end_time:
    
    stellar_II.evolve_model(model_time)
    channel_II.copy()
    
    print('Evolution at', model_time.value_in(units.Gyr), 'Gyr.' )
    model_time += step
    
stellar_II.stop()

# %%

# get the initial masses of the accreting stars
accreted_masses = initial_masses[mask[0]] | units.MSun
# make a particle set
rejuvenated_population = Particles(mass = accreted_masses)

# Evolve seperately the stars in the rejuvenated population in order to 
# define different metallicities for each one
for i in range(len(rejuvenated_population)):
    evolve_single_star(rejuvenated_population, i,
                       altered_metallicities[i], end_time, step)

# %%

HR_diagramme(old_population, rejuvenated_population, 20, save_to = directory_path)

# %%
