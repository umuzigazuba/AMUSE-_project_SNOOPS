#%%

from amuse.units import units
from amuse.lab import Particles, nbody_system, new_kroupa_mass_distribution
from amuse.community.seba.interface import SeBa
from amuse.ic.kingmodel import new_king_model

import numpy as np

from plotters import plot_snapshot_and_HR

# %%

def stellar_evolution(stars, metallicity, end_time, seed):
    '''
    Description:
        Evolve a star particle set to a desired age 
        
    Inputs:
        stars (Object): AMUSE particle set for a star system

        metallicity (float): Metallicity of the star system
        
        end_time (units.quantity): Total length of the evolution
        
        seed (int): Randomness of the function
    
    Returns:
        stars (Object): AMUSE particle set for the evolved star system
    '''
    np.random.seed(seed)

    stellar_evolution_code = SeBa()
    stellar_evolution_code.parameters.metallicity = metallicity
    stellar_evolution_code.particles.add_particles(stars)
    
    stellar_channel = stellar_evolution_code.particles.new_channel_to(stars)
    stellar_channel.copy()
    
    model_time = 0 | units.Myr
    dt = end_time / 50

    while(model_time < end_time):

        model_time += dt

        stellar_evolution_code.evolve_model(model_time)
        stellar_channel.copy()

    stellar_evolution_code.stop()
    return stars

#%%

def make_globular_cluster(star_count, radius, metallicity, age, seed):
    '''
    Description:
        Generate a globular cluster particle set

    Inputs:
        star_count (Int): Number of stars in the cluster
        
        radius (units.quantity): Core radius of the cluster
        
        metallicity (float): Metallicity of the stars
        
        age (units.quantity): Stellar evolution's timescale 

        seed (int): Randomness of the function
    
    Returns:
        evolved_cluster (object): AMUSE particle set for the globular cluster
    '''

    np.random.seed(seed)

    mass_kroupa = new_kroupa_mass_distribution(star_count,
                                                mass_min = 0.2 | units.MSun, 
                                                mass_max = 7 | units.MSun)
    stars = Particles(mass = mass_kroupa)

    evolved_stars = stellar_evolution(stars, metallicity, age, seed)
    converter = nbody_system.nbody_to_si(evolved_stars.mass.sum(), radius)


    evolved_cluster = new_king_model(star_count, W0 = 7, convert_nbody = converter)
    evolved_cluster.scale_to_standard(converter)
    
    evolved_cluster.mass = evolved_stars.mass
    evolved_cluster.luminosity = evolved_stars.luminosity
    evolved_cluster.temperature = evolved_stars.temperature
    evolved_cluster.metallicity = metallicity

    plot_snapshot_and_HR(evolved_cluster)
    
    return evolved_cluster

# %%
