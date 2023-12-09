##This file contains functions that are hard coded for the purpose of colliding a 
#globular cluster with a molecular cloud with our specific parameters 

from molecular_cloud_initialization import *
from plotters import *
from cluster_cloud_initialization import *

import os
import numpy as np
import matplotlib.pyplot as plt
from amuse.community.fi.interface import Fi
from amuse.lab import Particles, nbody_system
from amuse.couple import bridge
from amuse.units import units
from amuse.community.bhtree.interface import Bhtree
from amuse.ext.sink import new_sink_particles


def bondi_radius(stellar_mass):
    sound_speed = 0.2 | units.kms
    R = 2 * units.constants.G * stellar_mass /(sound_speed **2)

    return R

def bondi_accretion_rate(rho,v,r):
    dM = units.constants.pi * (r**2) * rho * v

    return dM

def detect_bounded_gas(star,particles,hardness):

    n = len(particles)

    if n == 0:
        return Particles()
    
    total_Ek=(0.5*particles.mass*(particles.vx**2+particles.vy**2+particles.vz**2)).sum()
    average_Ek=total_Ek/particles.mass.sum()
    limitE=hardness*average_Ek

    a=np.argsort(particles.x.number)
    binaries = Particles()

    for i in range(n):
        r2=(star.x-particles.x[a[i]])**2+ \
            (star.y-particles.y[a[i]])**2+ \
            (star.z-particles.z[a[i]])**2 
        v2=(star.vx-particles.vx[a[i]])**2+ \
            (star.vy-particles.vy[a[i]])**2+ \
            (star.vz-particles.vz[a[i]])**2 
        r=r2**0.5
        #Specific binding energy in units energy per mass
        eb=abs(units.constants.G*(particles.mass[a[i]]+star.mass)/r-0.5*v2)
        if eb > limitE:
            binary=particles[[a[i]]].copy()
            binary.hardness=eb/average_Ek
            binaries.add_particle(binary)
    
    return binaries

def free_fall_time(star, particles, time_step,bondi_radius):

    n = len(particles)
    
    if n == 0:
        return Particles()
    
    a = np.argsort(particles.x.number)
    binaries = Particles()
    
    rho = particles.mass.sum()/(4/3*np.pi*bondi_radius**3)

    for i in range(n):
        
        distance = np.sqrt((star.x-particles.x[a[i]])**2 + \
                            (star.y-particles.y[a[i]])**2 + \
                            (star.z-particles.z[a[i]])**2)
        
        tff = 0.5*np.pi*distance**(3/2)/np.sqrt(2*units.constants.G*(star.mass + particles.mass[a[i]]))

        if tff < time_step:
            binary = particles[[a[i]]].copy()
            binaries.add_particle(binary)

    return binaries


def accrete_mass(sinks, hydro_particles,time_step):
    # For each sink, find the hydro particles that are located within the sink radius
    particles_within_sink_radius = sinks.select_too_close(hydro_particles)

    for idx in range(len(sinks)):
        # Select the ones that are gravitationally bound to the sink
        bounded_particles = detect_bounded_gas(sinks[idx], particles_within_sink_radius[idx], hardness = 0.1)
        #bounded_particles = free_fall_time(sinks[idx], bounded_particles, time_step)
        if len(bounded_particles) != 0:
            # Update the mass of the sink
            sinks[idx].name = "Accreted star"
            sinks[idx].mass += bounded_particles.mass.sum()
            # Remove the accreted particles from the particle cloud
            hydro_particles.remove_particles(bounded_particles)


def make_cluster_with_vinit(velocity,position,random_seed):
    star_cluster = make_globular_cluster(star_count = 200,
                                        imf = "kroupa", 
                                        radius = 4 | units.pc,
                                        metallicity = 0.002, 
                                        age = 10 | units.Gyr, 
                                        seed = random_seed)
    star_cluster.name = "Unchanged star"
    
    print("Most massive star in cluster is ", max(star_cluster.mass.in_(units.MSun)))

    star_cluster.position +=  (-1.0, 0, 0) * (position | units.pc)
    star_cluster.velocity += (1.0, 0, 0) * (velocity| units.kms)

    converter_cluster = nbody_system.nbody_to_si(star_cluster.mass.sum(), 
                                    star_cluster.position.sum())
 
    return star_cluster,converter_cluster


#The following functions are hard coded for the convenience of running the
#simulation in batch

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


def AMUSE_bridge_initialization(star_cluster,converter_cluster,init_cloud,init_cloud_converter):
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



def let_them_collide_and_save(directory_path,t_end,dt,sinks,gravhydrobridge,\
                     channel,particles_cloud,gravity_code,hydro_cloud):
    t_end = t_end | units.Myr
    model_time = 0 | units.Myr
    dt = dt | units.Myr
    x_lim = int(abs(min(sinks.position.x).value_in(units.pc)))*1.2
    y_lim = int(abs(min(sinks.position.y).value_in(units.pc)))*2
    L=int(max(particles_cloud.x.value_in(units.pc)))*1.5
    N = 300

    sinks_mass_snapshots = []
    star_position = []
    cloud_density_cubes =[]

    density_map = plot_hydro(model_time,hydro_cloud,L,L,N)

    current_velocity = sum(sinks.vx.value_in(units.kms))/len(sinks.vx)
    print("Colliding with cluster velocity", current_velocity)

    while model_time < t_end:
        
        # define the accreting radius of the sinks particle based on its Bondi radius
        # IMPORTANT: the mass changes after each accretion event
        sinks.sink_radius = bondi_radius(sinks.mass).in_(units.pc)

        print("Largest sink radius", max(sinks.sink_radius).in_(units.pc))

        print("Pre accretion cluster mass", sinks.mass.sum().in_(units.MSun))
        
        model_time += dt
        model_time = model_time.round(1)
        # evolve the gravity and hydro codes through our bridge
        gravhydrobridge.evolve_model(model_time)


        # update channels (copy over from the codes.particles to the particle sets)
        channel["to_sinks"].copy()
        channel["to_cloud"].copy()


        print("Sinks in progress at", model_time.value_in(units.Myr), " Myr.")
        # add the acretted mass to the sinks's total mass

        accrete_mass(sinks,particles_cloud,dt)

        # update channels (copy the information from the particle set to the gravity code)
        channel["from_sinks"].copy()
        channel["from_cloud"].copy()

        # save necessary diagnostics of each step
        rho,_,_,_ = make_3Dmap(hydro_cloud,20,20)
        rho = rho.value_in(units.amu / units.cm**3)
        sinks_mass_snapshots.append(sinks.mass.value_in(units.MSun))
        star_position.append(sinks.position.value_in(units.pc))
        cloud_density_cubes.append(rho)

        if model_time+dt >= t_end:
            _, xgrid, ygrid, zgrid = make_3Dmap(hydro_cloud,20,20)

        print("Post accretion cluster mass", sinks.mass.sum().in_(units.MSun))
        print(len(particles_cloud.mass), "number of cloud particles now")

        density_plots_path = os.path.join(directory_path,"density_snapshots/")
        plot_cloud_and_star_cluster(model_time, hydro_cloud, sinks, x_lim, y_lim, N,density_map,saveto=density_plots_path)

    animation_path = os.path.join(directory_path, f"collision animation at{current_velocity}.html")
    fig = animate_collision_3D(star_position,cloud_density_cubes,xgrid,ygrid,zgrid)
    fig.write_html(animation_path)

    mass_difference = sinks_mass_snapshots[-1] - sinks_mass_snapshots[0]

    # Save the list to a text file
    with open(os.path.join(directory_path, \
                           f"Sink mass_{current_velocity}.txt"), 'w') as file:
        for sublist in sinks_mass_snapshots:
            line = ' '.join(map(str, sublist))  # Convert sublist to a space-separated string
            file.write(line + '\n')

    print("Mass snapshots saved.")

    mask = np.where(mass_difference > 1e-15)

    plt.plot(np.arange(0, t_end.value_in(units.Myr), dt.value_in(units.Myr)), \
             np.array(sinks_mass_snapshots)[:, mask[0]])
    plt.xlabel("time [Myr]")
    plt.ylabel("mass [Msun]")
    plt.savefig(os.path.join(directory_path, f"Sink mass with accretion_{current_velocity}.png"))
    plt.show()
    plt.close()

    relative_mass = np.array(mass_difference)[mask[0]]/np.array(sinks_mass_snapshots[0])[mask[0]]
    plt.hist(relative_mass*100, bins  = 30)
    plt.xlabel("Relative mass difference [%]")
    plt.savefig(os.path.join(directory_path, f"Accreted mass hist_{current_velocity}.png"))
    plt.close()

    final_cluster = sinks.copy()

    gravity_code.stop()
    hydro_cloud.stop()
    gravhydrobridge.stop()

    return final_cluster

