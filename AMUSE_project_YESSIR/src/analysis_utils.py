from amuse.units import units
from amuse.community.seba.interface import SeBa

import numpy as np


def updated_metallicity(M_star, M_material, Z_star = 0.002, Z_material = 0.02):
    '''
    Description:

        This function computes the metallicity of a star that has accreted some 
        material of mass M_material. The metallicity here is viewed as
        the mass fraction of elements hevier than Helium (Z = 1-X-Y).
    
    Inputs:
        
        M_star (float): The initial mass of the star
        
        M_material (float): The mass of the accreted material
        
        Z_star (float): The initial metallicity of the star (default value is
        the metallicity of Population II)
        
        Z_material (float): The metallicity of the accreted material (default 
        value is the metallicity of Population I)
        
    Outputs:
        
        Z_new (float): the new metallicity of the star, post accretion
    '''
    
    star_contribtion = (1-Z_star)*M_star
    material_contribution = (1-Z_material)*M_material
    M_new = M_star+M_material
    
    Z_new = 1 - np.divide(star_contribtion+material_contribution, M_new)
    # in cases of no accretion, machine error introduce a difference of ~1E-18
    Z_new = np.round(Z_new, 7)
    
    return Z_new


def evolve_single_star(stars, indx, metallicity, end_time, step):
    '''
    Description:
    
        This function evolves each star in a particle set seperately. The
        quantities of the evovled star is copied overe to the particle set
    
    Inputs:
        
        stars (amuse.datamodel.particles.Particles): The particle set
        containing the stars we want to evolve
        
        indx (int): the index f the particular star we want to evolve
        
        metallicity (float): the particular's star metallicity
        
        end_time (quantity): The duration of the simulation
        
        step (quantity): The time step of the simulation

    Outputs:

        None
    '''
    
    stellar = SeBa()
    # ALWAYS ADD THE METALLICITY FIRST AND THEN THE PARTICLES
    stellar.parameters.metallicity = metallicity
    stellar.particles.add_particle(stars[indx])
    channels = {"to_stars": stellar.particles.new_channel_to(stars), 
                "to_stellar": stars.new_channel_to(stellar.particles)}
    
    model_time = 0.0 | units.Gyr

    while model_time < end_time:
        
        stellar.evolve_model(model_time)
        channels["to_stars"].copy()
        
        #print('Evolution at', model_time.value_in(units.Gyr), 'Gyr' )
        model_time += step
        
    stellar.stop()


