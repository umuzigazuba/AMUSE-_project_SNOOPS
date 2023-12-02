#%%
import os
os.chdir('../src')

from molecular_cloud_initialization import *
from plotters import *
from cluster_cloud_initialization import *
from time import time
from datetime import datetime
import socket


import numpy as np
import matplotlib.pyplot as plt
from amuse.community.fi.interface import Fi
from amuse.lab import Particles, nbody_system
from amuse.couple import bridge
from amuse.units import units
from amuse.community.bhtree.interface import Bhtree
from amuse.ext.sink import new_sink_particles
#%%


tot_cloud_mass = 4/3 *units.constants.pi * (15 | units.pc)**3 * ( 2.3 | units.amu * 10 / (units.cm**3))
print(tot_cloud_mass.value_in(units.MSun))

#Assmuing a pure molecular hydrogen cloud, with typical density around 85 molecules per cm^-3, calculate the approximate cloud mass based on
#cloud size. 50 pc is selected for a small GC of only 100 stars. 

# %%
# initialise and evolve the MC particle set
init_cloud, init_cloud_converter  = make_molecular_cloud(N_cloud = 400_000,
                                                         M_cloud = 20_000 | units.MSun,
                                                         R_cloud = 15 | units.pc,
                                                         seed = 1312)

init_cloud, density_map = evolve_molecular_cloud(init_cloud, 
                                                    init_cloud_converter, 
                                                    t_end = 2 | units.Myr, 
                                                    dt = 0.2 | units.Myr, 
                                                    seed = 1312)

print("Mean density of the moelcular cloud", np.mean(init_cloud.density))
#%%
def hydro_code(Code, dt, converter, particles, seed):
    '''
    This function contains the parameters we want to initialise the 
    hydro code with. (hard Coded)
    '''
    
    np.random.seed(seed)

    hydro = Code(converter)
    hydro.parameters.use_hydro_flag = True # Hydrodynamics flag. True means:
                            # SPH hydro included, False means: gravity only.
    hydro.parameters.gamma = 1 # gas polytropic index (1.6666667)
                        # (default value:1.6666667). In this case-> Ideal Gas   
    hydro.parameters.timestep = dt
    hydro.parameters.eps_is_h_flag = True # Default value
    hydro.parameters.radiation_flag = False # turns off radiatiative cooling/heat.
    hydro.parameters.isothermal_flag = True  # Isothermal flag. True means:
                                            # isothermal gas (requires integrate_entropy_flag == False)
    hydro.parameters.integrate_entropy_flag = False #True means: integrate
                                                     # entropy, else: internal energy. 
    hydro.gas_particles.add_particles(particles) # add the particles
   
    return hydro    

def bondi_radius(stellar_mass):
    sound_speed = 0.2 | units.kms
    R = 2 * units.constants.G * stellar_mass /(sound_speed **2)

    return R

def bondi_accretion_rate(rho,v,r):
    dM = units.constants.pi * (r**2) * rho * v

    return dM


def make_cluster_with_vinit(velocity):
    star_cluster = make_globular_cluster(star_count = 200,
                                        imf = "kroupa", 
                                        radius = 4 | units.pc,
                                        metallicity = 0.002, 
                                        age = 10 | units.Gyr, 
                                        seed = 2804)

    star_cluster.position +=  (-1.0, 0, 0) * (30 | units.pc)
    star_cluster.velocity += (1.0, 0, 0) * (velocity| units.kms)

    converter_cluster = nbody_system.nbody_to_si(star_cluster.mass.sum(), 
                                    star_cluster.position.sum())
 
    return star_cluster,converter_cluster





def AMUSE_bridge_initialization(star_cluster,converter_cluster):
    #initiate the gravity code with sink particles
    gravity_code = Bhtree(converter_cluster)
    sinks = new_sink_particles(star_cluster)

    gravity_code.particles.add_particles(sinks)

    # #start the hydro code for the gas
    converter_cloud = init_cloud_converter
    particles_cloud = init_cloud.copy()
    hydro_cloud = hydro_code(Code = Fi, dt = 0.1 | units.Myr,
                            converter = converter_cloud,
                            particles = particles_cloud,
                            seed = 1312)


    channel = {"to_sinks": gravity_code.particles.new_channel_to(sinks),
                "from_sinks": sinks.new_channel_to(gravity_code.particles),
                "to_cloud": hydro_cloud.gas_particles.new_channel_to(particles_cloud),
                "from_cloud": particles_cloud.new_channel_to(hydro_cloud.gas_particles)}


    gravhydrobridge = bridge.Bridge(use_threading = False)
    gravhydrobridge.add_system(gravity_code, (hydro_cloud,) )
    gravhydrobridge.add_system(hydro_cloud, (gravity_code,) )
    gravhydrobridge.timestep = 0.1 | units.Myr

    return gravhydrobridge,sinks,channel,particles_cloud,gravity_code,hydro_cloud


def let_them_collide(directory_path,t_end,dt,sinks,gravhydrobridge,channel,particles_cloud,gravity_code,hydro_cloud):
    t_end = t_end | units.Myr
    model_time = 0 | units.Myr
    dt = dt | units.Myr

    sinks_mass_snapshots = []
    current_velocity = sum(sinks.vx.value_in(units.kms))/len(sinks.vx)
    print("Colliding with cluster velocity", current_velocity)

    while model_time < t_end:
        
        # define the accreting radius of the sinks particle based on its Bondi radius
        # IMPORTANT: the mass changes after each accretion event
        sinks.sink_radius = bondi_radius(sinks.mass).in_(units.pc)

        print("Largest sink radius", max(sinks.sink_radius).in_(units.pc))

        print("Pre accretion cluster mass", sinks.mass.sum())
        
        model_time += dt
        model_time = model_time.round(1)
        # evolve the gravity and hydro codes through our bridge
        gravhydrobridge.evolve_model(model_time)


        # update channels (copy over from the codes.particles to the particle sets)
        channel["to_sinks"].copy()
        channel["to_cloud"].copy()

        print("Sinks in progress at", model_time.value_in(units.Myr), " Myr.")
        # add the acretted mass to the sinks's total mass
        sinks.accrete(particles_cloud)

        # update channels (copy the information from the particle set to the gravity code)
        channel["from_sinks"].copy()
        channel["from_cloud"].copy()

        # save the total mass of each step
        sinks_mass_snapshots.append(sinks.mass.value_in(units.MSun))

        print("Post accretion cluster mass", sinks.mass.sum().in_(units.MSun))
        print(len(particles_cloud.mass), "number of cloud particles now")


        plt.scatter(particles_cloud.x.value_in(units.pc), particles_cloud.y.value_in(units.pc), s = 1)
        plt.scatter(sinks.x.value_in(units.pc), sinks.y.value_in(units.pc), c = 'red', s = 5)
        plt.title("Molecular cloud at time = " + model_time.as_string_in(units.Myr))
        plt.xlabel("x [pc]")
        plt.ylabel("y [pc]")
        plt.savefig(os.path.join(directory_path, f"at time_{model_time.value_in(units.Myr)}.png"))
        plt.close()


    mass_difference = sinks_mass_snapshots[-1] - sinks_mass_snapshots[0]

    # Save the list to a text file
    with open(os.path.join(directory_path, \
                           f"Sink mass_{current_velocity}.txt"), 'w') as file:
        for sublist in sinks_mass_snapshots:
            line = ' '.join(map(str, sublist))  # Convert sublist to a space-separated string
            file.write(line + '\n')

    print("Mass snapshots saved.")

    plt.plot(sinks_mass_snapshots)
    plt.xlabel("time [Myr]")
    plt.ylabel("mass [Msun]")
    plt.show()
    plt.savefig(os.path.join(directory_path, f"Sink mass_{current_velocity}.png"))
    plt.close()

    plt.hist(mass_difference, bins  = 30)
    plt.savefig(os.path.join(directory_path, f"Accreted mass hist_{current_velocity}.png"))
    plt.close()

    gravity_code.stop()
    hydro_cloud.stop()
    gravhydrobridge.stop()

#%%


velocity = np.linspace(15,70,10,dtype=int)

for i in range (len(velocity)):
    v = velocity[i]
    print("Starting with cluster velocity",v)
    star_cluster,converter_cluster = make_cluster_with_vinit(v)
    gravhydrobridge,sinks,channel,particles_cloud,gravity_code,\
        hydro_cloud = AMUSE_bridge_initialization(star_cluster,converter_cluster)
    
    path = '../results/'
    variable_name = v
    directory_path = os.path.join(path, str(variable_name))
    os.makedirs(directory_path, exist_ok=True)

    if __name__ == "__main__":
        start_time = time()
        date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
        print("Collision started on {}".format(socket.gethostname()),f"at time {format(date_time)}")
        
        let_them_collide(directory_path,2,0.1,sinks,gravhydrobridge,channel,particles_cloud,gravity_code,hydro_cloud)

        print("Collision finished (running time: {0:.1f}s)".format(time() - start_time))
# %%
