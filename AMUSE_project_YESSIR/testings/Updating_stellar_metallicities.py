

import numpy as np
import matplotlib.pyplot as plt

from amuse.units import units
from amuse.community.seba.interface import SeBa
from amuse.datamodel import Particles
#%%
collision_velocity = 20
directory_path = '../results/final_with_stellar_evolution/20 kms/'
mass_file = 'Sink mass_19.977764918756957.txt'
# import the mass snapshots file (in Solar Mass)
sinks_mass_snapshots = np.loadtxt(directory_path + mass_file)
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
    # in cases of no accretion, machine error introduce a difference of ~1E-18
    Z_new = np.round(Z_new, 7)
    
    return Z_new



new_metallicities = updated_metallicity(M_star= initial_masses,
                                        M_material=accreted_mass)


#%%

# mask the elements that show change in their metallicity
mask = np.where(new_metallicities > 0.002)

altered_metallicities = new_metallicities[mask[0]]


#%%
# Z HISTOGRAM
alpha = 1
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_facecolor('whitesmoke')
ax.grid(alpha=alpha/2)

ax.hist(altered_metallicities, bins  = 30,
        color = '#0c2577', label =  'altered Z')
ax.axvline(x=0.002, linestyle = '--', color = 'k', label = 'initial Z')


fig.supylabel('#', size = 'x-large')
fig.supxlabel('Z', size = 'x-large')
fig.suptitle('Collision with velocity '+ str(collision_velocity) +' kms',
             size = 'x-large')
ax.legend(fontsize = 14)
plt.tight_layout()
plt.savefig(directory_path+'metallicity_hist.png')

#%%
# make the mask for the non-accreting stars
mask_non_accretion = np.where(new_metallicities == 0.002)
# get the initial masses of the non accreting stars
non_accreted_masses = initial_masses[mask_non_accretion[0]] | units.MSun
# make a particle set of the old population
old_population = Particles(mass = non_accreted_masses)

# Stellar evolution of old population
 

stellar_II = SeBa()
# ALWAYS ADD THE METALLICITY FIRST AND THEN THE PARTICLES
stellar_II.parameters.metallicity = 0.002
stellar_II.particles.add_particles(old_population)
channel_II = stellar_II.particles.new_channel_to(old_population)
channel_II.copy()


# evovle for the age of the cluster + the time spend in the collision with the
# molecular cloud
end_time = (10 |units.Gyr) +(1.8 |units.Myr)
model_time = 0.0 |units.Gyr
step = 2 |units.Gyr


while model_time < end_time:
    
    stellar_II.evolve_model(model_time)
    channel_II.copy()
    
    print('Evolution at', model_time.value_in(units.Gyr), 'Gyr' )
    model_time += step
    
stellar_II.stop()

#%%
# get the initial masses of the accreting stars
accreted_masses = initial_masses[mask[0]] | units.MSun
# make a particle set
rejuvenated_population = Particles(mass =accreted_masses)
#%%
def evolve_single_star(stars, indx, metallicity, end_time, step):
    '''
    This function evolves each star in a particle set seperately. The
    quantities of the evovled star is copied overe to the particle set
    
    Inputs:
        
        stars (amuse.datamodel.particles.Particles): The particle set
        containing the stars we want to evolve
        
        indx (int): the index f the particular star we want to evolve
        
        metallicity (float): the particular's star metallicity
        
        end_time (quantity): The duration of the simulation
        
        step (quantity): The time step of the simulation
    '''
    
    stellar = SeBa()
    # ALWAYS ADD THE METALLICITY FIRST AND THEN THE PARTICLES
    stellar.parameters.metallicity = metallicity
    stellar.particles.add_particle(stars[indx])
    channels = {"to_stars": stellar.particles.new_channel_to(stars), 
                "to_stellar": stars.new_channel_to(stellar.particles)}
    

    print('The metallicity of the current star is Z = ',stellar.parameters.metallicity)
    model_time = 0.0 |units.Gyr

    while model_time < end_time:
        
        stellar.evolve_model(model_time)
        channels["to_stars"].copy()
        
        #print('Evolution at', model_time.value_in(units.Gyr), 'Gyr' )
        model_time += step
        
    stellar.stop()



#%%

# Evolve seperately the stars in the rejuvenated population in order to 
# define different metallicities for each one
for i in range(len(rejuvenated_population)):
    
    evolve_single_star(rejuvenated_population, i,
                       altered_metallicities[i], end_time, step)



#%%

# HR DIAGRAMME
fig, ax = plt.subplots( figsize=(9,6))
ax.set_facecolor('whitesmoke')
ax.grid(alpha=alpha/2)

ax.scatter(old_population.temperature.value_in(units.K),
                old_population.luminosity.value_in(units.LSun), 
                c='k', s= 12,
                label = 'Stars without accretion')


ax.scatter(rejuvenated_population.temperature.value_in(units.K),
                  rejuvenated_population.luminosity.value_in(units.LSun), 
                  c='red', s=6,
                  label = 'Stars with accretion')

ax.set_xlim(8E+3, 3.15E+3)
ax.set_ylim(1E-4, 5E0)

plt.axvspan(3150, 3500, facecolor='sandybrown', alpha=0.7, zorder = 0)
plt.axvspan(3500, 5000, facecolor='navajowhite', alpha=0.7, zorder = 0)
plt.axvspan(5000, 6000, facecolor='khaki', alpha=0.7, zorder = 0)
plt.axvspan(6000, 7500, facecolor='lightyellow', alpha=0.7, zorder = 0)
plt.axvspan(7500, 8000, facecolor='lightcyan', alpha=0.7, zorder = 0)
# MS TRACK, ZOOMED IN
#ax.set_xlim(6.5E+3, 3.15E+3)
#ax.set_ylim(1E-3, 5E0)
########################################

x1, x2, y1, y2 = 3.55E+3, 4E+3, 5E-3, 5E-2  # subregion of the original image
axins = ax.inset_axes(
    [.15, 0.3, .2, .4],
    xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

axins.set_facecolor('whitesmoke')

axins.scatter(old_population.temperature.value_in(units.K),
                old_population.luminosity.value_in(units.LSun), 
                c='k', s= 12,
                label = 'old population')

axins.scatter(rejuvenated_population.temperature.value_in(units.K),
                  rejuvenated_population.luminosity.value_in(units.LSun), 
                  c='red', s=12,
                  label = 'rejuvenated population')


axins.set_xlim(4E+3, 3.55E+3)
axins.get_xaxis().set_visible(False)
axins.set_ylim(5E-3, 5E-2)
axins.get_yaxis().set_visible(False)
axins.loglog()




##################################################
ax.indicate_inset_zoom(axins, edgecolor="black")
ax.loglog()
ax.legend(fontsize = 14)
ax.set_xlabel("T [K]", size = 'x-large')
ax.set_ylabel(r"L [L$_{\odot}$]", size = 'x-large')
fig.suptitle('Collision with velocity '+ str(collision_velocity) +' kms',
             size = 'x-large')
plt.tight_layout()
plt.savefig(directory_path+'HR_spctrl.png')

