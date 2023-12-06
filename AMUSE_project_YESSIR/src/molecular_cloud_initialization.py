# %%
from amuse.lab import nbody_system

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.community.fi.interface import Fi

from amuse.units import units
import numpy as np
from matplotlib import pyplot as plt

from plotters import plot_hydro

# %%

def make_molecular_cloud(N_cloud, M_cloud, R_cloud, seed):
    '''
    Description:
        Generating a molecular cloud
        
    Inputs:
        N_cloud (Int): Number of particles in the molecular cloud
        
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

def evolve_molecular_cloud(particles_cloud, converter_cloud, t_end, dt, seed):
    '''
    Description:
        Evolve a existing molecular cloud to a certain age

    Inputs:
        particles_cloud (object): AMUSE particle set for the molecular cloud 

        converter_cloud (object): AMUSE generic unit converter

        t_end (units.quantity): Total length of the evolution

        dt (units.quantity): Time step of the evolution

        seed (int): Randomness of the function 

    Return:
        particles_cloud (object): AMUSE particle set for the evolved molecular cloud

        density_map (matplotlib.image.AxesImage): AxesImage object containing the plotted log density map of the evolved molecular cloud

    '''

    np.random.seed(seed)

    hydro_cloud = Fi(converter_cloud)
    hydro_cloud.gas_particles.add_particles(particles_cloud)

    hydro_cloud.parameters.use_hydro_flag = True # Hydrodynamics flag. True means: SPH hydro included, False means: gravity only.
    hydro_cloud.parameters.radiation_flag = False # Radiation flag. True means: 
                                                #radiation (i.e. radiative cooling/heating) is included. 
    hydro_cloud.parameters.timestep = dt # timestep for system (default value:4.70451599238e+13 s)

    hydro_cloud.parameters.gamma = 1 # gas polytropic index (1.6666667) (default value:1.6666667)
    hydro_cloud.parameters.isothermal_flag = True # Isothermal flag. True means: isothermal gas (requires integrate_entropy_flag == False). (default value:False)
    hydro_cloud.parameters.integrate_entropy_flag = False # Integrate-entropy flag. True means: integrate entropy, else: internal energy. (default value:True)
       
    channel = {"hydro_to_part": hydro_cloud.gas_particles.new_channel_to(particles_cloud),
            "part_to_hydro": particles_cloud.new_channel_to(hydro_cloud.gas_particles)}

    L = int(max(particles_cloud.x.value_in(units.pc)))*1.5  # x and y lim of plot. 
    N = 1000 # amount of grid points

    model_time = 0 | units.Myr

    density_map = plot_hydro(model_time, hydro_cloud, L, L, N)
    
    print("ready for evolution")
    while model_time < t_end:

        model_time += dt
        model_time = model_time.round(2)

        hydro_cloud.evolve_model(model_time)
        print("Time", model_time.in_(units.Myr))
        channel["hydro_to_part"].copy()

    density_map = plot_hydro(model_time, hydro_cloud, L, L, N)

    hydro_cloud.stop()

    print(f"Average mass of a SPH particle {particles_cloud.mass.sum().value_in(units.MSun)/len(particles_cloud.mass)}.")

    return particles_cloud, density_map

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

