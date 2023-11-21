#%%
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath("AMUSE_project_YESSIR"))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from src.molecular_cloud_initialization import make_molecular_cloud, \
                                            evolve_molecular_cloud, make_map

#%%
import numpy as np
import matplotlib.pyplot as plt
from amuse.community.fi.interface import Fi
from amuse.lab import Particles, nbody_system
from amuse.couple import bridge
from amuse.units import units
from amuse.community.bhtree.interface import Bhtree
#%%

particles_cloud, converter_cloud  = make_molecular_cloud(N_cloud = 50_000,
                                                         M_cloud = 1_000 | units.MSun,
                                                         R_cloud = 20 | units.parsec,
                                                         seed = 1312)
particles_cloud = evolve_molecular_cloud(particles_cloud, 
                                         converter_cloud, 
                                         2 | units.Myr, 
                                         0.2 | units.Myr, 
                                         1 | units.RSun,
                                         1312)


#%%
collision_bodies = Particles(1)
collision_bodies[0].name = "star"
collision_bodies[0].mass = 25 |units.MSun
collision_bodies[0].radius = 8 | units.RSun
collision_bodies[0].position = (-1.0,0,0) * (45 | units.pc)
collision_bodies[0].velocity = (1.0,0,0) * (60 | units.kms)
converter_star=nbody_system.nbody_to_si(collision_bodies[0].mass.sum(), 
                                   collision_bodies[0].position.length())

gravity_code = Bhtree(converter_star)
gravity_code.particles.add_particles(collision_bodies)
# ch_gravity2star = gravity_code.particles.new_channel_to(collision_bodies)

#%%
resolution = 1 | units.RSun
timestep = 0.2 | units.Myr

np.random.seed(1312)

hydro_code = Fi(converter_cloud)

hydro_code.parameters.use_hydro_flag = True # Hydrodynamics flag. True means: 
                                             #SPH hydro included, False means: gravity only.
hydro_code.parameters.radiation_flag = False # Radiation flag. True means: 
                                                #radiation (i.e. radiative cooling/heating) is included. 

hydro_code.parameters.gamma = 1 # gas polytropic index (1.6666667) (default value:1.6666667)
hydro_code.parameters.isothermal_flag = True # Isothermal flag. True means: isothermal gas (requires integrate_entropy_flag == False). (default value:False)
hydro_code.parameters.integrate_entropy_flag = False # Integrate-entropy flag. True means: integrate entropy, else: internal energy. (default value:True)
hydro_code.parameters.timestep = 0.5*timestep # timestep for system (default value:4.70451599238e+13 s)

hydro_code.parameters.eps_is_h_flag = False  #Eps-is-h flag. True means: set gas particles gravitational epsilon to h (SPH smoothing length). (default value:True)
                                            # h_smooth is constant
hydro_code.parameters.gas_epsilon = resolution # The gas gravitational smoothing epsilon.
hydro_code.parameters.sph_h_const = resolution # SPH smoothing length if constant
particles_cloud.h_smooth= resolution

hydro_code.gas_particles.add_particles(particles_cloud)



#%%
star = collision_bodies[collision_bodies.name=='star']
hydro_code.dm_particles.add_particles(star)
collision_bodies.add_particles(particles_cloud)
channel = {"from star": star.new_channel_to(gravity_code.particles),
            "to_star": gravity_code.particles.new_channel_to(star),
            "from_cloud": particles_cloud.new_channel_to(hydro_code.particles),
            "to_cloud":hydro_code.particles.new_channel_to(particles_cloud),
            "hydro_from_star": star.new_channel_to(hydro_code.dm_particles),
            "hydro_to_star": hydro_code.dm_particles.new_channel_to(star)}


# collision_bodies.add_particles(particles_cloud)
#%%
gravhydrobridge = bridge.Bridge(use_threading=False)
gravhydrobridge.add_system(gravity_code, (hydro_code,) )
gravhydrobridge.add_system(hydro_code, (gravity_code,) )
gravhydrobridge.timestep = 0.2|units.Myr
#%%

def plot_cloud_star(time, hydro_code, star_particle, L,N):
    fig = plt.figure(figsize = (9, 5))
 
    rho = make_map(hydro_code, L = L, N = N)
    cax = plt.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), extent=[-L,L,-L,L]) # , vmin = 0, vmax = 5
    plt.scatter(star_particle.x.value_in(units.pc), star_particle.y.value_in(units.pc), c = 'red')
    cbar = fig.colorbar(cax)
    cbar.set_label('log density [$amu/cm^3$]', labelpad = 5)

    plt.title("Molecular cloud at time = " + time.as_string_in(units.Myr))
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()

def zoomed_in_star(time, hydro_code, star_particle, L, x_center,y_center,N):
    fig = plt.figure(figsize = (9, 5))
 
    rho = make_map(hydro_code, L = L, N = N)
    cax = plt.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), extent=[-L,L,-L,L]) # , vmin = 0, vmax = 5
    
    cbar = fig.colorbar(cax)
    cbar.set_label('log density [$amu/cm^3$]', labelpad = 5)

    offset = L/10 + 5

    plt.axis([x_center-offset, x_center+offset, y_center-5, y_center+5])
    plt.scatter(star_particle.x.value_in(units.pc), star_particle.y.value_in(units.pc), c = 'red')
    plt.title("Zoomed in molecular cloud at time = " + time.as_string_in(units.Myr))
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()
#%%

model_time = 0 | units.Myr
dt = timestep
t_end = 1.4 | units.Myr
L = int(max(particles_cloud.x.value_in(units.pc))) 
L = L + L

while model_time < t_end:

    model_time += dt


    gravhydrobridge.evolve_model(model_time)
    channel["to_star"].copy()
    channel["to_cloud"].copy()
    channel["hydro_to_star"].copy()
    print('Collision in progress at ',model_time.value_in(units.Myr), 'Myr')
    print("star position", star.position.value_in(units.pc), "pc")
    plot_cloud_star(model_time, hydro_code, collision_bodies[0], L,1000)

    star_x= star.position[0][0].value_in(units.pc)
    star_y= star.position[0][1].value_in(units.pc)

    zoomed_in_star(model_time, hydro_code, collision_bodies[0], L, star_x,star_y,1000)

gravity_code.stop()
hydro_code.stop()



# %%
print("star position", star.position[0][1].value_in(units.pc))
print("star mass", star.mass)
# %%
plot_cloud_star(t_end, hydro_code, collision_bodies[0], 5, 1000)
# %%