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

particles_cloud, converter_cloud  = make_molecular_cloud(N_cloud = 1_00,
                                                         M_cloud = 1_000 | units.MSun,
                                                         R_cloud = 10 | units.parsec,
                                                         seed = 1312)

particles_cloud, colorbar = evolve_molecular_cloud(particles_cloud, 
                                         converter_cloud, 
                                         t_end = 2 | units.Myr, 
                                         dt = 0.2 | units.Myr, 
                                         seed = 1312)


#%%
collision_bodies = Particles(1)
collision_bodies[0].name = "star"
collision_bodies[0].mass = 150 |units.MSun
collision_bodies[0].radius = 35 | units.RSun
collision_bodies[0].position = (-1.0,0,0) * (20 | units.pc)
collision_bodies[0].velocity = (1.0,0,0) * (20 | units.kms)
converter_star=nbody_system.nbody_to_si(collision_bodies[0].mass.sum(), 
                                   collision_bodies[0].position.length())

gravity_code = Bhtree(converter_star)
gravity_code.particles.add_particles(collision_bodies)

#%%
timestep = 0.1 | units.Myr

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

def plot_cloud_star(time, hydro_code, star_particle, L, x_center, y_center, N, colorbar):
    rho = make_map(hydro_code, L=L, N=N)
    
    fig, (ax_full, ax_zoom) = plt.subplots(nrows = 2, ncols = 1, figsize=(10.5, 7))
 
    im_full = ax_full.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), extent=[-L, L, -L, L])
    ax_full.scatter(star_particle.x.value_in(units.pc), star_particle.y.value_in(units.pc), c='red')

    ax_full.set_title("Molecular cloud at time = " + time.as_string_in(units.Myr))
    ax_full.set_xlabel("x [pc]")
    ax_full.set_ylabel("y [pc]")

    im_zoom = ax_zoom.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), extent=[-L, L, -L, L])
    ax_zoom.scatter(star_particle.x.value_in(units.pc), star_particle.y.value_in(units.pc), c='red')

    offset = L/10 + 5

    ax_zoom.axis([x_center - offset, x_center + offset, y_center - 5, y_center + 5])
    ax_zoom.set_title("Zoomed in molecular cloud at time = " + time.as_string_in(units.Myr))
    ax_zoom.set_xlabel("x [pc]")
    ax_zoom.set_ylabel("y [pc]")
    
    cbar_ax = fig.add_axes([0.75, 0.1, 0.02, 0.85])
    cbar = fig.colorbar(colorbar, cax = cbar_ax)
    cbar.set_label('log density [$amu/cm^3$]', labelpad = 5)

    plt.tight_layout()
    plt.show()

#%%

model_time = 0 | units.Myr
dt = timestep
t_end = 2 | units.Myr
L = int(max(particles_cloud.x.value_in(units.pc))) 
L = L + L

star_x_position = []
star_y_position = []
star_z_position = []

while model_time < t_end:

    model_time += dt
    model_time = model_time.round(2)

    gravhydrobridge.evolve_model(model_time)
    channel["to_star"].copy()
    channel["to_cloud"].copy()
    channel["hydro_to_star"].copy()
    
    print('Collision in progress at ',model_time.value_in(units.Myr), 'Myr')

    star_x= star.position[0][0].value_in(units.pc)
    star_y= star.position[0][1].value_in(units.pc)
    star_z= star.position[0][2].value_in(units.pc)

    star_x_position.append(star_x)
    star_y_position.append(star_y)
    star_z_position.append(star_z)

    plot_cloud_star(model_time, hydro_code, collision_bodies[0], L, star_x, star_y, 500, colorbar)

gravity_code.stop()
hydro_code.stop()

# %%
plt.scatter(star_x_position, star_y_position)
plt.title(f"(x, y)-position of star over {t_end.value_in(units.Myr)} Myr.")
plt.xlabel("x [pc]")
plt.ylabel("y [pc]")
plt.show()

plt.scatter(star_x_position, star_z_position)
plt.title(f"(x, z)-position of star over {t_end.value_in(units.Myr)} Myr.")
plt.xlabel("x [pc]")
plt.ylabel("z [pc]")
plt.show()

# %%
# print("star position", star.position[0][1].value_in(units.pc))
# print("star mass", star.mass)
# # %%
# plot_cloud_star(t_end, hydro_code, collision_bodies[0], 5, 1000)
# # %%
