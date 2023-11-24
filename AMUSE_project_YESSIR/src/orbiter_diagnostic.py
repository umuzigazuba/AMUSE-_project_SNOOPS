import numpy as np
import matplotlib.pyplot as plt
from amuse.community.fi.interface import Fi
from amuse.lab import Particles, nbody_system
from amuse.couple import bridge
from amuse.units import units
from amuse.community.bhtree.interface import Bhtree
import os
os.chdir('/home/kostas/AMUSE/AMUSE_project_YESSIR/AMUSE_project_YESSIR/src')

from molecular_cloud_initialization import *
from cluster_cloud_initialization import *

#%%

def hydro_code(Code, dt, converter, particles, seed):
    '''
    This function contains the parameters we want to initialise the 
    hydro code with. (hard Coded)
    '''
    
    np.random.seed(seed)

    hydro = Fi(converter)
    hydro.parameters.use_hydro_flag = True # Hydrodynamics flag. True means:
                            # SPH hydro included, False means: gravity only.
    hydro.parameters.gamma = 1 # gas polytropic index (1.6666667)
                        # (default value:1.6666667). In this case-> Ideal Gas   
    hydro.parameters.timestep = dt
    hydro.parameters.eps_is_h_flag = True # Default value
    hydro.parameters.radiation_flag = False # turns off radiatiative cooling/heat.
    hydro.parameters.isothermal_flag = True  # Isothermal flag. True means:
                # isothermal gas (requiresintegrate_entropy_flag == False)
    hydro.parameters.integrate_entropy_flag = False #True means: integrate
                                          # entropy, else: internal energy. 
    hydro.gas_particles.add_particles(particles) # add the particles
   
    return hydro 

# initialise the gas particle set
particles_cloud, converter_cloud  = make_molecular_cloud(N_cloud = 1_000,
                                                         M_cloud = 100 | units.MSun,
                                                         R_cloud = 20 | units.parsec,
                                                         seed = 1312)

#initialise the star particle

ALL_bodies = Particles(1)
ALL_bodies[0].name = "star"
ALL_bodies[0].mass = 100 |units.MSun
ALL_bodies[0].radius = 10 | units.RSun
ALL_bodies[0].position = (1.0,0,0) * (45 | units.pc)
ALL_bodies[0].velocity = (0,1.0,0) * (20 | units.kms)
converter_star=nbody_system.nbody_to_si(ALL_bodies[0].mass.sum(), 
                                   ALL_bodies[0].position.length())




#%%
#start the hydro code for the gas
hydro_cloud = hydro_code(Code = Fi, dt = 0.2 |units.Myr,
                         converter = converter_cloud,
                         particles = particles_cloud,
                         seed = 1312)

# start the gravity code for the star
gravity_code = Bhtree(converter_star)
gravity_code.particles.add_particles(ALL_bodies)
#%%

star = ALL_bodies[ALL_bodies.name=='star']

ALL_bodies.add_particles(particles_cloud)

#%%


channel = {"from_star": ALL_bodies.new_channel_to(gravity_code.particles),
            "to_star": gravity_code.particles.new_channel_to(ALL_bodies),
            "from_cloud": ALL_bodies.new_channel_to(hydro_cloud.particles),
            "to_cloud":hydro_cloud.particles.new_channel_to(ALL_bodies)}
#            "hydro_from_star": ALL_bodies.new_channel_to(hydro_cloud.dm_particles),
#            "hydro_to_star": hydro_cloud.dm_particles.new_channel_to(ALL_bodies)}

#%%

gravhydrobridge = bridge.Bridge(use_threading=False)
gravhydrobridge.add_system(gravity_code, (hydro_cloud,) )
gravhydrobridge.add_system(hydro_cloud, (gravity_code,) )
gravhydrobridge.timestep = 0.2|units.Myr

#%%
##############################################
# WRONG MASSES AND VELOCITY FOR ORBIT. BAD DEFINITION 
# OF CLOUD'S CENTER OF MASS (IT ALSO INCLUED THE MASS/POS OF STAR
# AND THE CENTER OF MASS IS DRIFTED TOWARDS THE STAR'S POS)

t_end = 250 |units.Myr
model_time = 0 |units.Myr
dt = 5 |units.Myr
X_star = [] | units.pc
Y_star = [] | units.pc
Z_star = [] | units.pc

X_cloud = [] | units.pc
Y_cloud = [] | units.pc
Z_cloud = [] | units.pc
while model_time < t_end:

    model_time += dt
    gravhydrobridge.evolve_model(model_time)

    channel["to_star"].copy()
    channel["to_cloud"].copy()
    # channel["hydro_to_star"].copy()
    # save the star's position
    X_star.append(ALL_bodies[0].x)
    Y_star.append(ALL_bodies[0].y)
    Z_star.append(ALL_bodies[0].z)
    # save the cloud's position
    x_cloud = hydro_cloud.particles.center_of_mass()[0]
    X_cloud.append(x_cloud)
    y_cloud = hydro_cloud.particles.center_of_mass()[1]
    Y_cloud.append(y_cloud)
    z_cloud = hydro_cloud.particles.center_of_mass()[2]
    Z_cloud.append(z_cloud)
    print('Orbiter in progress at',model_time.value_in(units.Myr), ' Myr')

#%%
gravity_code.stop()
hydro_cloud.stop()
gravhydrobridge.stop()
#%%


from amuse.plot import plot
plot(X_star, Y_star, lw=1)
plot(X_cloud, Y_cloud)
#pyplot.gca().set_aspect("equal", adjustable="box")
#pyplot.show()
