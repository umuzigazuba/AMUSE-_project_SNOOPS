#%%
from amuse.units import units
from amuse.community.seba.interface import SeBa
import numpy as np
import matplotlib.pyplot as plt
from amuse.ext.masc import new_star_cluster

# %%
def stellar_evolution(bodies, time, random_seed):
    np.random.seed(random_seed)
    stellar_evolution_code = SeBa()
    stellar_evolution_code.particles.add_particles(bodies)
    stellar_evolution_code.parameters.metallicity = bodies[0].metallicity
    
    stellar_channel = stellar_evolution_code.particles.new_channel_to(bodies)
    stellar_channel.copy()
    

    end_time = time
    model_time = 0 | units.Myr
    dt = end_time/50
    print(dt.in_(units.Myr))
    while(model_time<end_time):
        np.random.seed(random_seed)
        model_time += dt
        stellar_evolution_code.evolve_model(model_time)
        stellar_channel.copy()
  
    stellar_evolution_code.stop()
    plot_snapshot(bodies)
    return bodies

#%%
def plot_snapshot(bodies):
    v = (bodies.vx**2 + bodies.vy**2 + bodies.vz**2).sqrt()
    s = bodies.mass.value_in(units.MSun)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3))
    ax1.scatter(bodies.temperature.value_in(units.K),
                bodies.luminosity.value_in(units.LSun), 
                c='k',
                s= 10)
    ax1.set_xlim(7E+3, 2.5E+3)
    #ax1.set_ylim(1.e+3, 1.e+7)
    ax1.loglog()
    ax1.set_xlabel("T [K]")
    ax1.set_ylabel("L [$L_\odot$]")
    ax2.scatter(bodies.x.value_in(units.pc), 
                bodies.y.value_in(units.pc), 
                c=v.value_in(units.kms), 
                s=s)
    plt.gca().set_aspect('equal', adjustable='box')
    ax2.set_xlabel("x [pc]")
    ax2.set_ylabel("y [pc]")
    #ax2.set_xlim(-5, 5)
    #ax2.set_ylim(-5, 5)
    plt.show()

#%%

def make_globular_cluster(star_count, imf, metallicity, age, random_seed):
    '''
    Description:
        
    Inputs:
        star_count (Int): Number of stars in the cluster
        
        imf (str): Either 'kroupa' or 'salpeter'
        
        metallicity (float): The stars' metallicity
        
        age (units.quantity): Stellar evolution's timescale 
    '''
    np.random.seed(random_seed)
    cluster = new_star_cluster(
        number_of_stars=star_count,
        initial_mass_function=imf,
        upper_mass_limit= 100 | units.MSun, 
        effective_radius= 50 | units.parsec, #assuming this is the overall radius of the cluster
        star_distribution='king',
        star_distribution_w0=7.0,
        star_metallicity=metallicity,
    )
    print("cluster generated")

    evolved_cluster = stellar_evolution(cluster, age, random_seed)

    plot_snapshot(cluster)
    
    return evolved_cluster

# %%

new_cluster = make_globular_cluster(1000,'kroupa',0.002,3 ,723476)

# %%





