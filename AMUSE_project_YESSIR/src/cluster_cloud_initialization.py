#%%
from amuse.units import units
from amuse.community.seba.interface import SeBa
import numpy as np
import matplotlib.pyplot as plt
from amuse.ext.masc import new_star_cluster

from plotters import plot_snapshot_and_HR

# %%
def stellar_evolution(bodies, time, seed):
    '''
    Description:
        Evolve an existing star particles system to a desired age 
        
    Inputs:
        bodies (Object): AMUSE particle set for a star system
        
        time (units.quantity): Evolution time scale 
        
        seed (int): Randomness of the function
    
    Returns:
        bodies (Object): AMUSE particle set for the evolved star system
    '''
    np.random.seed(seed)

    stellar_evolution_code = SeBa()
    stellar_evolution_code.particles.add_particles(bodies)
    stellar_evolution_code.parameters.metallicity = bodies[0].metallicity
    
    stellar_channel = stellar_evolution_code.particles.new_channel_to(bodies)
    stellar_channel.copy()
    

    end_time = time
    model_time = 0 | units.Myr
    dt = end_time / 50

    while(model_time < end_time):

        model_time += dt
        stellar_evolution_code.evolve_model(model_time)
        stellar_channel.copy()
  
    stellar_evolution_code.stop()
    return bodies

#%%

def make_globular_cluster(star_count, imf, radius, metallicity, age, seed):
    '''
    Description:
        Generate a globular cluster with specific number of stars,
         inital mass function, metallicity, and age.

    Inputs:
        star_count (Int): Number of stars in the cluster
        
        imf (str): Initial mass function. Either 'kroupa' or 'salpeter'

        radius (units.quantity): The effective radius of the cluster
        
        metallicity (float): The stars' metallicity
        
        age (units.quantity): Stellar evolution's timescale 

        seed (int): Randomness of the function
    
    Returns:
        evolved_cluster (object): AMUSE particle system resembling the desired globular cluster
    '''

    np.random.seed(seed)

    cluster = new_star_cluster(
        number_of_stars = star_count,
        initial_mass_function = imf,
        upper_mass_limit = 100 | units.MSun, 
        effective_radius = radius, #assuming this is the overall radius of the cluster
        star_distribution = 'king',
        star_distribution_w0 = 7.0,
        star_metallicity = metallicity,
    )

    print("cluster generated")

    evolved_cluster = stellar_evolution(cluster, age, seed)

    plot_snapshot_and_HR(evolved_cluster)
    
    return evolved_cluster

# %%

#new_cluster = make_globular_cluster(1000,'kroupa',0.002,3|units.Gyr,723476)
# IMPORTANT: The age has a large effect on the amount of evolved stars present (e.g. 3 vs 10 Gyr) 
# QUESTION: evolve until the average age difference between populations in globular clusters
# %%




