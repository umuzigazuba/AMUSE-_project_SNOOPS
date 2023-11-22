# %%
from amuse.lab import nbody_system

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.community.fi.interface import Fi

from amuse.units import units
import numpy as np
from matplotlib import pyplot as plt

# %%

# plot molecular cloud density function (smooth)
def make_map(hydro, L, N):
    '''
    Description:

    Inputs:
        hydro (object): AMUSE Fi hydrodynamic integrator
        L (int): Axis length 
        N (int): Number of grid points

    Return: 
        rho ()
    '''

    x = np.linspace(-L, L, N + 1)
    y = np.linspace(-L, L, N + 1)
    xv, yv = np.meshgrid(x, y)

    x = xv.flatten() | units.pc
    y = yv.flatten() | units.pc
    z = 0 | units.pc
    vx = 0 | units.kms
    vy = 0 | units.kms
    vz = 0 | units.kms

    rho = hydro.get_hydro_state_at_point(x, y, z, vx, vy, vz)[0]
    rho = rho.reshape((N + 1, N + 1))
    
    return rho
    
def plot_hydro(time, hydro, L, N):

    fig = plt.figure(figsize = (9, 5))
    
    rho = make_map(hydro, L = L, N = N)
    cax = plt.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), extent=[-L, L, -L, L])
    cbar = fig.colorbar(cax)
    cbar.set_label('log density [$amu/cm^3$]', labelpad = 5)
        
    plt.title("Molecular cloud at time = " + time.as_string_in(units.Myr))
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()

    return cax

def plot_cloud_particles(time, particles_cloud):

    fig = plt.figure(figsize = (9, 5))

    plt.scatter(particles_cloud.x.value_in(units.pc), particles_cloud.y.value_in(units.pc), s = 1)
    plt.title("Molecular cloud at time = " + time.as_string_in(units.Myr))
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()

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
        resolution (units.quantity): Gas particle smoothing length 
        seed (int): Randomness of the function 

    Return:
        particles_cloud (object): AMUSE particle set for the evolved molecular cloud

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

    initial_colorbar = plot_hydro(model_time, hydro_cloud, L, N)
    print("ready for evolution")

    while model_time < t_end:

        model_time += dt
        model_time = model_time.round(1)

        hydro_cloud.evolve_model(model_time)
        print("Time", model_time.in_(units.Myr))
        channel["hydro_to_part"].copy()

    final_colorbar = plot_hydro(model_time, hydro_cloud, L, N)

    hydro_cloud.stop()

    print(f"Average mass of a SPH particle {particles_cloud.mass.sum().value_in(units.MSun)/len(particles_cloud.mass)}.")

    return particles_cloud, final_colorbar

# %%

# particles_cloud, converter_cloud = make_molecular_cloud(100_000, 10_000_000 | units.MSun, 200 | units.parsec, 1312)
# particles_cloud = evolve_molecular_cloud(particles_cloud, converter_cloud, 2 | units.Myr, 0.2 | units.Myr, 1312)
