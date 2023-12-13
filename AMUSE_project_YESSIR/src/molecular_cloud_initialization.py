# %%

from amuse.lab import nbody_system
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.community.fi.interface import Fi
from amuse.units import units

import numpy as np

from plotters import plot_hydro

# %%

def make_molecular_cloud(N_cloud, M_cloud, R_cloud, seed):
    '''
    Description:
        Generate a molecular cloud particle set
        
    Inputs:
        N_cloud (Int): Target number of particles in the molecular cloud
        
        M_cloud (units.quantity): Total mass of the cloud
        
        R_cloud (units.quantity): Radius of the cloud
        
        seed (int): Randomness of the function
    
    Returns:
        particles_cloud (object): AMUSE particle set for the molecular cloud

        converter_cloud (object): AMUSE generic unit converter for the cloud 
    '''

    converter_cloud = nbody_system.nbody_to_si(M_cloud, R_cloud)

    # creates a smooth spherical cloud with random velocities as in Bonnell et al. (2003)
    particles_cloud = molecular_cloud(targetN = N_cloud, 
                                        convert_nbody = converter_cloud,
                                        base_grid = body_centered_grid_unit_cube,
                                        seed = seed).result 
    
    return particles_cloud, converter_cloud

# %%

def evolve_molecular_cloud(particles_cloud, converter_cloud, end_time, dt, seed):
    '''
    Description:
        Evolve an existing molecular cloud to a certain age

    Inputs:
        particles_cloud (object): AMUSE particle set for the molecular cloud 

        converter_cloud (object): AMUSE generic unit converter

        end_time (units.quantity): Total length of the evolution

        dt (units.quantity): Time step of the evolution

        seed (int): Randomness of the function 

    Return:
        particles_cloud (object): AMUSE particle set for the evolved molecular cloud

        density_map (matplotlib.image.AxesImage): AxesImage object containing the plotted log density map of the evolved molecular cloud

    '''

    np.random.seed(seed)

    hydro = Fi(converter_cloud)

    hydro.parameters.use_hydro_flag = True # SPH hydrodynamics is included alongside self-gravity
    hydro.parameters.radiation_flag = False # Radiation flag  is included

    hydro.parameters.timestep = dt # Timestep used by the code

    hydro.parameters.gamma = 1 # Gas polytropic index
    hydro.parameters.isothermal_flag = True # Isothermal gas
    hydro.parameters.integrate_entropy_flag = False # Integrate entropy

    hydro.gas_particles.add_particles(particles_cloud)
       
    channel = {"hydro_to_part": hydro.gas_particles.new_channel_to(particles_cloud),
               "part_to_hydro": particles_cloud.new_channel_to(hydro.gas_particles)}


    L = int(max(particles_cloud.x.value_in(units.pc)))*1.5  # X- and y-limit of the density plot
    N = 1000 # Amount of grid points used in the density plot

    model_time = 0 | units.Myr

    density_map = plot_hydro(model_time, hydro, L, L, N)
    
    print("Ready for evolution")
    while model_time < end_time:

        model_time += dt
        model_time = model_time.round(1)
        print("Time", model_time.in_(units.Myr))

        hydro.evolve_model(model_time)
        channel["hydro_to_part"].copy()

    density_map = plot_hydro(model_time, hydro, L, L, N)

    hydro.stop()

    print(f"Average mass of a SPH particle {particles_cloud.mass.sum().value_in(units.MSun) / len(particles_cloud.mass)}.")

    return particles_cloud, density_map

# %%
