# This file contains functions that are hard coded for the purpose of colliding a 
# globular cluster with a molecular cloud with our specific parameters 

# %%
 
from cluster_initialization import make_globular_cluster
from plotters import plot_cloud_and_star_cluster, plot_evolution_mass_accretion, plot_relative_mass

from amuse.units import units
from amuse.lab import Particles, nbody_system
from amuse.ext.sink import new_sink_particles
from amuse.community.bhtree.interface import Bhtree
from amuse.community.sse.interface import SSE
from amuse.community.fi.interface import Fi
from amuse.couple import bridge

import os
import numpy as np
import matplotlib.pyplot as plt

# %%

def bondi_radius(stellar_mass):
    '''
    Description:

        Calculate the bondi radius of a star

    Input:

        stellar_mass (units.quantity): Mass of the star

    Outputs:

        bondi_radius (units.quantity): Bondi radius of the star
    '''

    sound_speed = 0.2 | units.kms # Sound speed in a typical molecular cloud
    bondi_radius = 2 * units.constants.G * stellar_mass / sound_speed**2

    return bondi_radius

# %%

def detect_bounded_gas(star, particles, hardness):
    '''
    Description:

        Return the particles that are gravitationally bound to a star

    Input:

        star (object): AMUSE sink particle 

        particles (object): AMUSE particle set containing the molecular cloud
        particles that lie within the sink radius of the star

        hardness (float): Harshness of the gravitationally bound criterion

    Outputs:

        bounded_particles (object): AMUSE particle set containing the molecular
        cloud particles that are gravitationally bound to the star
    '''

    N = len(particles)

    if N == 0:
        return Particles()
    # compute the energy limit to detect agravitationally bound binary
    total_kinetic_energy = (0.5 * particles.mass * (particles.vx**2 + particles.vy**2 + particles.vz**2)).sum()
    average_kinetic_energy = total_kinetic_energy / particles.mass.sum()
    limit_energy = hardness * average_kinetic_energy

    sorted_idx = np.argsort(particles.x.number)
    bounded_particles = Particles()
    # evaluate each pair of star - i-th gas particle
    for i in range(N):

        r2 = (star.x - particles.x[sorted_idx[i]])**2 + \
             (star.y - particles.y[sorted_idx[i]])**2 + \
             (star.z - particles.z[sorted_idx[i]])**2 
        
        v2 = (star.vx - particles.vx[sorted_idx[i]])**2 + \
             (star.vy - particles.vy[sorted_idx[i]])**2 + \
             (star.vz - particles.vz[sorted_idx[i]])**2 
        
        r = r2**0.5

        #Specific binding energy in units energy per mass
        binding_enegry = abs(units.constants.G * (particles.mass[sorted_idx[i]] + star.mass) / r - 0.5*v2)

        if binding_enegry > limit_energy:
            valid_particle = particles[[sorted_idx[i]]].copy()
            bounded_particles.add_particle(valid_particle)
    
    return bounded_particles

# %%

def free_fall_time(star, particles, bounded_particles, time_step):
    '''
    Description:

        Return the particles that are gravitationally bound to a star 
        and how the amount of mass that is accreted per particle by the star
        within one time step

        The amount of mass depends on the free-fall time

    Inputs:

        star (object): AMUSE sink particle 

        particles (object): AMUSE particle set containing the molecular cloud
        particles that lie within the sink radius of the star

        bounded_particles (object): AMUSE particle set containing the
        gravitationally bound molecular cloud particles

        time_step (units.quantity): Time interval within which mass is accreted

    Outputs:

        bounded_particles (object): AMUSE particle set containing the
        gravitationally bound molecular cloud particles

        accreted_mass (numpy.ndarray): The amount of mass accreted from each
        gravitationally bound particle
    '''

    N = len(bounded_particles)
    accreted_mass = []
    
    if N == 0:
        return Particles(), accreted_mass
    
    sorted_idx = np.argsort(bounded_particles.x.number)
    accreted_particles = Particles()
    
    background_density = particles.mass.sum() / (4 * np.pi * star.sink_radius**3 / 3)

    for i in range(N):
        
        distance = ((star.x - bounded_particles.x[sorted_idx[i]])**2 + \
                    (star.y - bounded_particles.y[sorted_idx[i]])**2 + \
                    (star.z - bounded_particles.z[sorted_idx[i]])**2)**0.5
        
        free_fall_radius = (2 * time_step * np.sqrt(2 * units.constants.G \
                              * (star.mass + bounded_particles.mass[sorted_idx[i]])) / np.pi)**(2/3)
        
        valid_particle = bounded_particles[[sorted_idx[i]]].copy()
        accreted_particles.add_particle(valid_particle)

        if free_fall_radius >= distance:
            accreted_mass.append(valid_particle.mass.value_in(units.MSun)[0])

        else:
            valid_particle_denisty = valid_particle.mass / (4 * np.pi * distance**3 / 3)
            mean_rho = (0.5 * background_density + 1.5 * valid_particle_denisty) / 2
            acquired_mass = (mean_rho * 4 * np.pi * free_fall_radius**3) / 3 
            accreted_mass.append(acquired_mass.value_in(units.MSun)[0])
    
    accreted_mass = accreted_mass | units.MSun
    return accreted_particles, accreted_mass

# %%

def accrete_mass(sinks, hydro_particles, time_step):
    '''
    Description:

        Determine the accretion of molecular cloud particles by a star cluster
        within a certain timestep

    Inputs:

        sinks (object): AMUSE sink particle set represeting the star cluster

        hydro_particles (object): AMUSE particle set represeting the molecular
        cloud particles

        time_step (units.quantity): Time interval within which mass is accreted

    Outputs:

        None
    '''

    # For each sink, find the hydro particles that are located within the sink radius
    adjacent_particles = sinks.select_too_close(hydro_particles)

    for idx in range(len(sinks)):

        # Select the particles that are gravitationally bound to the sink
        bounded_particles = detect_bounded_gas(sinks[idx], adjacent_particles[idx], hardness = 0.1)
        # Calculate the amount of mass that is accreted from each particle
        bounded_particles, accreted_mass = free_fall_time(sinks[idx],
                                                adjacent_particles[idx],
                                                bounded_particles, time_step)
        
        if len(bounded_particles) != 0:
            # Update the mass of the sink
            sinks[idx].name = "Accreted star"
            sinks[idx].mass += np.sum(accreted_mass)

            #Update the mass of the molecular cloud particles
            bounded_particles.mass -= accreted_mass

            for i in range(len(bounded_particles)):

                if bounded_particles[i].mass.value_in(units.MSun) <= 1e-15:
                    hydro_particles.remove_particle(bounded_particles[i])

                else:
                    for particles in hydro_particles:   

                        if particles.key == bounded_particles.key[i]:
                            particles.mass = bounded_particles.mass[i]
                    
# %%

def make_cluster_with_posvel(position, velocity, W0, random_seed,
                             number_of_stars = 200):
    '''
    Description:

        Initialize a globular cluster, its position and velocity 

    Inputs:

        position (float): Position of the star cluster

        velocity (float): Velocity of the star cluster

        W0 (float): King model parameter

        random_seed (int): Randomness of the cluster initialization

        number_of_stars (int): Number of stars in the star cluster

    Outputs:

        star_cluster (object): AMUSE particel set for the star cluster

        converter_cluster (object): AMUSE generic unit converter for the star
        cluster
    '''

    star_cluster = make_globular_cluster(star_count = number_of_stars, 
                                         radius = 4 | units.pc,
                                         metallicity = 0.002, 
                                         age = 10 | units.Gyr, 
                                         W0 = W0,
                                         seed = random_seed)
    
    star_cluster.name = "Unchanged star"

    print(f"Least massive star in the cluster has a mass {min(star_cluster.mass.in_(units.MSun))}.")  
    print(f"Most massive star in the cluster has a mass {max(star_cluster.mass.in_(units.MSun))}.")

    star_cluster.position += (-1.0, 0, 0) * (position | units.pc)
    star_cluster.velocity += (1.0, 0, 0) * (velocity | units.kms)

    converter_cluster = nbody_system.nbody_to_si(star_cluster.mass.sum(), 
                                                 star_cluster.position.sum())
 
    return star_cluster, converter_cluster

# %%

def hydrodynamics_code(code, time_step, particles_cloud, converter_cloud, seed):
    '''
    Description:

        Initialize the hydrodynamics code for a molecular cloud
    
    Inputs:

        code (object): AMUSE hydrodynamics code

        time_step (units.quantity): Time step of the hydrodynamics code

        particles_cloud (object): AMUSE particle set for the molecular cloud

        converter_cloud (object): AMUSE generic unit converter for the
        molecular cloud

        seed (int): Randomness of the function

    Outputs:

        hydro (object): AMUSE hydrodynamics integrator for the molecular cloud
    '''
    
    np.random.seed(seed)

    hydro = code(converter_cloud)

    hydro.parameters.use_hydro_flag = True # SPH hydrodynamics is included alongside self-gravity
    hydro.parameters.radiation_flag = False # Radiation flag  is included

    hydro.parameters.timestep = time_step # Timestep used by the code

    hydro.parameters.gamma = 1 # Gas polytropic index
    hydro.parameters.isothermal_flag = True # Isothermal gas
    hydro.parameters.integrate_entropy_flag = False # Integrate entropy

    hydro.gas_particles.add_particles(particles_cloud)
   
    return hydro    
##############################################################################
# The following functions are hard coded for the convenience of running the
# simulation in batch
##############################################################################


def code_bridge_channel_initaization(time_step, star_cluster, converter_cluster,
                                     particles_cloud, converter_cloud, seed):
    '''
    Description:

        Initialize the codes, channels and bridge for the collsion of a star
        cluster and molecular cloud
    
    Inputs:

        time_step (units.quantity): Time step of the hydrodynamics code

        star_cluster (object): AMUSE particel set for the star cluster

        converter_cluster (object): AMUSE generic unit converter for the star
        cluster

        particles_cloud (object): AMUSE particle set for the molecular cloud

        converter_cloud (object): AMUSE generic unit converter for the
        molecular cloud

        seed (int): Randomness of the function
    
    Outputs:

        sinks (object): AMUSE sink particle set for the star cluster

        gravity_code (object): AMUSE gravity integrator for the star cluster

        stellar_evolution_code (object): AMUSE stellar evolution integrator for
        the star cluster

        hydro_code (object): AMUSE hydrodynamics integrator for the molecular
        cloud

        channels (dictionary): Channels between the particle sets and
        code.particles

        gravity_hydro_bridge (object): AMUSE bridge integrator between the
        gravity and hydrodynamics codes
    '''

    # Create a sink particle set with the same properties as the star cluster
    sinks = new_sink_particles(star_cluster)

    # Initiate the gravity code for the star cluster
    gravity_code = Bhtree(converter_cluster)
    gravity_code.particles.add_particles(sinks)

    # Initiate the stellar evolution code for the star cluster
    stellar_evolution_code = SSE()
    stellar_evolution_code.parameters.metallicity = sinks[0].metallicity
    stellar_evolution_code.particles.add_particles(sinks)

    # Initate the hydrodynamics code for the molecular cloud
    hydro_code = hydrodynamics_code(code = Fi, 
                                    time_step = time_step,
                                    converter_cloud = converter_cloud,
                                    particles_cloud = particles_cloud,
                                    seed = seed)

    # Initate the channels between the codes and particle sets
    channels = {"gravity_to_sinks": gravity_code.particles.new_channel_to(sinks),
                "gravity_from_sinks": sinks.new_channel_to(gravity_code.particles, attributes=["mass", "radius"], target_names=["mass", "radius"]),
                "stellar_evolution_to_sinks": stellar_evolution_code.particles.new_channel_to(sinks),
                "stellar_evolution_from_sinks": sinks.new_channel_to(stellar_evolution_code.particles),
                "hydro_to_cloud": hydro_code.gas_particles.new_channel_to(particles_cloud),
                "hydro_from_cloud": particles_cloud.new_channel_to(hydro_code.gas_particles)}

    # Update the star cluster properties to the stellar evolution properties  
    channels["stellar_evolution_to_sinks"].copy()

    # Initate the bridge between the gravity and hydrodynamics codes
    gravity_hydro_bridge = bridge.Bridge(use_threading = False)
    gravity_hydro_bridge.add_system(gravity_code, (hydro_code,))
    gravity_hydro_bridge.add_system(hydro_code, (gravity_code,))
    gravity_hydro_bridge.timestep = 2 * time_step

    return sinks, gravity_code, stellar_evolution_code, hydro_code, channels, gravity_hydro_bridge 

# %%

def cluster_cloud_collision(end_time, time_step, sinks, particles_cloud, gravity_code, stellar_evolution_code, \
                            hydro_code, channels, gravity_hydro_bridge, directory_path, density_map):
    '''
    Description:

        Perform the collision of a star cluster and molecular cloud and save the results 

    Inputs:

        end_time (units.quantity): Duration of the collision

        time_step (units.quantity): Time step of the collision

        sinks (object): AMUSE sink particle set for the star cluster

        particles_cloud (object): AMUSE particle set for the molecular cloud

        gravity_code (object): AMUSE gravity integrator for the star cluster

        stellar_evolution_code (object): AMUSE stellar evolution integrator for the star cluster

        hydro_code (object): AMUSE hydrodynamics integrator for the molecular cloud

        channels (dictionary): Channels between the particle sets and code.particles

        gravity_hydro_bridge (object): AMUSE bridge integrator between the gravity and hydrodynamics codes

        directory_path (str): Path to the folder where the results should be saved

        density_map (matplotlib.image.AxesImage): AxesImage object for the log density map of the molecular cloud before collision

    Outputs:

        final_cluster (object): AMUSE sink particle set for the star cluster after the collision
    '''

    sinks_mass_evolution = []

    current_velocity = np.mean(sinks.vx.value_in(units.kms))
    print(f"Colliding with cluster velocity {current_velocity}")

    model_time = 0 | units.Myr

    # Define the parameters for the density plot 
    x_limit = int(abs(min(sinks.position.x).value_in(units.pc)))*1.2
    y_limit = int(abs(min(sinks.position.y).value_in(units.pc)))*2
    N = 300

    density_plots_path = os.path.join(directory_path,"density_snapshots/")
    plot_cloud_and_star_cluster(time = model_time, hydro = hydro_code, sinks = sinks, \
                                x_limit = x_limit, y_limit = y_limit, N = N, density_map_MC = density_map, \
                                save_to = density_plots_path)

    while model_time < end_time:

        print(f"Collision time = {model_time.in_(units.Myr)}.")
              
        model_time += time_step
        model_time = model_time.round(1)

        stellar_evolution_code.evolve_model(model_time)

        # Update the sink particles and gravity_code.particles through the channels 
        channels["stellar_evolution_to_sinks"].copy()
        channels["gravity_from_sinks"].copy()
        
        sinks.sink_radius = bondi_radius(sinks.mass).in_(units.pc)

        print(f"Largest sink radius {max(sinks.sink_radius).in_(units.pc)}.")

        # Evolve the gravity and hydrodynamics codes through the bridge
        gravity_hydro_bridge.evolve_model(model_time)

        # Update the particle sets through the channels
        channels["gravity_to_sinks"].copy()
        channels["hydro_to_cloud"].copy()

        print(f"Pre accretion cluster mass {sinks.mass.sum().in_(units.MSun)}.")

        # Add the accreted mass to the sinks's total mass
        accrete_mass(sinks, particles_cloud, time_step)

        # Update codes.particles through the channels
        channels["stellar_evolution_from_sinks"].copy()
        channels["gravity_from_sinks"].copy()
        channels["hydro_from_cloud"].copy()

        print(f"Post accretion cluster mass {sinks.mass.sum().in_(units.MSun)}. \n")

        sinks_mass_evolution.append(sinks.mass.value_in(units.MSun))

        plot_cloud_and_star_cluster(time = model_time, hydro = hydro_code, sinks = sinks, \
                                    x_limit = x_limit, y_limit = y_limit, N = N, density_map_MC = density_map, \
                                    save_to = density_plots_path)

    # Save the stellar mass evolution to a text file
    with open(os.path.join(directory_path, f"Sink_mass_{current_velocity}.txt"), 'w') as file:

        for sublist in sinks_mass_evolution:
            line = ' '.join(map(str, sublist))  # Convert sublist to a space-separated string
            file.write(line + '\n')

    print("Mass snapshots saved.")

    plot_evolution_mass_accretion(sinks_mass_evolution = sinks_mass_evolution, 
                                  end_time = end_time, 
                                  time_step = time_step, 
                                  velocity = current_velocity, 
                                  save_to = directory_path)

    plot_relative_mass(sinks_mass_evolution = sinks_mass_evolution, 
                       velocity = current_velocity, 
                       save_to = directory_path)

    final_cluster = sinks.copy()

    gravity_code.stop()
    stellar_evolution_code.stop()
    hydro_code.stop()
    gravity_hydro_bridge.stop()

    return final_cluster

# %%