import os
os.chdir('./src')

import numpy as np
import matplotlib.pyplot as plt


from amuse.units import nbody_system
from amuse.units import units
from amuse.community.seba.interface import SeBa
#%%

file_path = '../results/testng_paths/45/Sink mass_44.97792473385582.txt'
# import the mass snapshots file (in Solar Mass)
sinks_mass_snapshots = np.loadtxt(file_path)
# get the masses of the Cluster's stars before and after the collision
initial_masses = sinks_mass_snapshots[0]
final_masses = sinks_mass_snapshots[-1]
# compute how much mass is acreeted from the molecular cloud
accreted_mass = np.subtract(final_masses, initial_masses)
#%%
def updated_metallicity(M_star, M_material, Z_star = 0.002, Z_material = 0.02):
    '''
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
    
    return Z_new



new_metallicities = updated_metallicity(M_star= initial_masses,
                                        M_material=accreted_mass)




