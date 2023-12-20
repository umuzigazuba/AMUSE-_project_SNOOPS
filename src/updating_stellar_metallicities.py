# %%

import numpy as np

from amuse.units import units
from amuse.community.seba.interface import SeBa
from amuse.datamodel import Particles

from analysis_utils import updated_metallicity, evolve_single_star
from plotters import metallicity_histogram, HR_diagramme

# %%

collision_velocity = 20
directory_path = '../results/final_with_stellar_evolution/20 kms/'
mass_file = 'Sink mass_19.977764918756957.txt'
# Import the mass snapshots file (in Solar Mass)
sinks_mass_snapshots = np.loadtxt(directory_path + mass_file)
# Get the masses of the Cluster's stars before and after the collision
initial_masses = sinks_mass_snapshots[0]
final_masses = sinks_mass_snapshots[-1]
# Compute how much mass is acreeted from the molecular cloud
accreted_mass = np.subtract(final_masses, initial_masses)
# Compute the new metallicities of all stars in the cluster
new_metallicities = updated_metallicity(M_star = initial_masses,
                                        M_material = accreted_mass)


# Mask the elements that show change in their metallicity
mask = np.where(new_metallicities > 0.002)
altered_metallicities = new_metallicities[mask[0]]


metallicity_histogram(altered_metallicities, N_bins = 30,
                      collision_velocity = collision_velocity)

# %%

# Make the mask for the non-accreting stars
mask_non_accretion = np.where(new_metallicities == 0.002)
# Get the initial masses of the non accreting stars
non_accreted_masses = initial_masses[mask_non_accretion[0]] | units.MSun
# Make a particle set of the old population
old_population = Particles(mass = non_accreted_masses)

# Stellar evolution of old population
stellar_II = SeBa()
# ALWAYS ADD THE METALLICITY FIRST AND THEN THE PARTICLES
stellar_II.parameters.metallicity = 0.002
stellar_II.particles.add_particles(old_population)
channel_II = stellar_II.particles.new_channel_to(old_population)
channel_II.copy()

# Evovle for the age of the cluster + the time spend in the collision with the
# molecular cloud
end_time = (10 | units.Gyr) + (1.8 | units.Myr)
model_time = 0.0 | units.Gyr
step = 2 | units.Gyr

while model_time < end_time:
    
    stellar_II.evolve_model(model_time)
    channel_II.copy()
    
    print('Evolution at', model_time.value_in(units.Gyr), 'Gyr' )
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

HR_diagramme(old_population, rejuvenated_population, 20)














