# %%

import os
from amuse.units import units

import numpy as np
from matplotlib import pyplot as plt


# %%

def plot_snapshot_and_HR(cluster, save_to = None):
    '''
    Description: 
        Plot the positions of stars in a star cluster at z = 0 
        and the Hertzsprung-Russell diagram of the cluster

    Inputs:
        cluster (object): AMUSE particle set for the star cluster

        save_to (str): Path to the folder where the image should be saved

    Return: 
        None
    '''

    velocity = (cluster.vx**2 + cluster.vy**2 + cluster.vz**2).sqrt()
    mass = cluster.mass

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9, 3))
    # HR diagramme
    ax1.scatter(cluster.temperature.value_in(units.K),
                cluster.luminosity.value_in(units.LSun), 
                c = 'k',
                s = 8)
    ax1.set_xlim(8.5E+3, 2.5E+3) 
    ax1.set_ylim(1.e-6, 1.e+3)
    ax1.loglog()
    ax1.set_xlabel("Temperature [K]")
    ax1.set_ylabel("Luminosity [$L_\odot$]")
    # x-y plane. colour based on velocities
    plot = ax2.scatter(cluster.x.value_in(units.pc), 
                        cluster.y.value_in(units.pc), 
                        c = velocity.value_in(units.kms), 
                        s = mass.value_in(units.MSun))
    plt.gca().set_aspect('equal', adjustable = 'box')
    ax2.set_xlabel("x [pc]")
    ax2.set_ylabel("y [pc]")

    colorbar = fig.colorbar(plot)
    colorbar.set_label('velocity [km/s]', fontsize = 9)

    if save_to is not None:
        os.makedirs(save_to, exist_ok = True)
        plt.savefig(os.path.join(save_to, f"Snapshot_and_HR_of_cluster.png"))
        plt.close()
    
    else:
        plt.show()
# %%

def make_map(hydro, x_limit, y_limit, N):
    '''
    Description: 
        For a molecular cloud initially centered around x, y, z = 0 pc, 
        calculate its density over a grid at z = 0 pc

    Inputs:
        hydro (object): AMUSE hydrodynamic integrator for the molecular cloud

        x_limit (int): Half the x-axis length in parsecs

        y_limit (int): Half the y-axis length in parsecs

        N (int): Number of grid points along one axis

    Return: 
        density (numpy.ndarray): Two-dimensional array for the density of the
        molecular cloud over the grid
    '''
    # create a grid on the x-y plane
    x = np.linspace(-x_limit, x_limit, N)
    y = np.linspace(-y_limit, y_limit, N)
    X, Y = np.meshgrid(x, y)

    x = X.flatten() | units.pc
    y = Y.flatten() | units.pc
    z = 0 | units.pc
    vx = 0 | units.kms
    vy = 0 | units.kms
    vz = 0 | units.kms
    # evaluate the cloud's density at the verices of the grid
    density = hydro.get_hydro_state_at_point(x, y, z, vx, vy, vz)[0]
    density = density.reshape((N, N))
    
    return density

# %%  

def plot_hydro(time, hydro, x_limit, y_limit, N, save_to = None):
    '''
    Description: 
        Plot the log density of a molecular cloud at a given time and z = 0 pc
        Used to track the pre-collision evolution of a molecular cloud 

    Inputs: 
        time (units.quantity): Time since the molecular cloud initialization

        hydro (object): AMUSE hydrodynamic integrator for the molecular cloud

        x_limit (int): Half the x-axis length in parsecs

        y_limit (int): Half the y-axis length in parsecs

        N (int): Number of grid points along one axis

        save_to (str): Path to the folder where the image should be saved

    Return:
        density_map (matplotlib.image.AxesImage): AxesImage object of the
        plotted log density map
    '''

    fig = plt.figure(figsize = (9, 5))
    
    density = make_map(hydro, x_limit = x_limit, y_limit = y_limit, N = N)
    density_map = plt.imshow(np.log10(density.value_in(units.amu/units.cm**3)), cmap = "plasma", extent = [-x_limit, x_limit, -y_limit, y_limit])
    
    color_bar = fig.colorbar(density_map)
    color_bar.set_label('log density [$amu/cm^3$]', labelpad = 5)
        
    plt.title(f"Molecular cloud at time = {time.value_in(units.Myr)} Myr and z = 0 pc")
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    
    if save_to is not None:
        os.makedirs(save_to, exist_ok = True)
        plt.savefig(os.path.join(save_to, f"Molecular_cloud_at_time_{time.value_in(units.Myr)}_Myr.png"))
        plt.close()
    
    else:
        plt.show()

    return density_map

# %%

def plot_cloud_and_star_cluster(time, hydro, sinks, x_limit, y_limit,
                                N, density_map_MC, save_to = None):
    '''
    Description: 
        Plot the log density of a molecular cloud and the position of a
        colliding star cluster at a given time Used to track the collision
        between a molecular cloud and a star cluster

    Inputs:
        time (units.quantity): Time since the start of the collision

        hydro (object): AMUSE hydrodynamic integrator containing the molecular
        cloud particles

        sinks (object): AMUSE sink particle set for the star cluster

        x_limit (int): Half the x-axis length in parsecs

        y_limit (int): Half the y-axis length in parsecs

        N (int): Number of grid points along one axis

        density_map_MC (matplotlib.image.AxesImage): AxesImage object for the
        log density map of the molecular cloud before collision

        save_to (str): Path to the folder where the image should be saved

    Return:
        None
    '''
    
    rho = make_map(hydro, x_limit = x_limit, y_limit = y_limit, N = N)

    colors = np.array(["black", "aliceblue"])
    colors_sink = np.array([0 if sink.name == "Unchanged star" else 1 \
                            for sink in sinks])

    fig = plt.figure(figsize = (10, 5))

    plt.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), cmap = "plasma",
               extent = [-x_limit, x_limit, -y_limit, y_limit])
    
    plt.scatter(sinks.position.x.value_in(units.pc),
                sinks.position.y.value_in(units.pc), c = colors[colors_sink],
                s = sinks.mass.value_in(units.MSun)*5)

    plt.title(f"Collision at time = {time.value_in(units.Myr)} Myr")
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.xlim([-x_limit, x_limit])
    plt.ylim([-y_limit, y_limit])

    colorbar_axis = fig.add_axes([0.85, 0.1, 0.02, 0.85])
    colorbar = plt.colorbar(density_map_MC, cax = colorbar_axis,
                            fraction = 0.046, pad = 0.04)
    colorbar.set_label('log density [$amu/cm^3$]', labelpad = 5)

    if save_to is not None:
        os.makedirs(save_to, exist_ok=True)
        plt.savefig(os.path.join(save_to, f"Molecular_cloud_and_cluster_at_time_{time.value_in(units.Myr)}_Myr.png"))
        plt.close()
    
    else:
        plt.show()

# %%

def plot_evolution_mass_accretion(sinks_mass_evolution, end_time, time_step,
                                  velocity, save_to = None):
    '''
    Description:
        Plot the amount of accreted mass as a function of time for each star
        in a cluster that accretes mass

    Inputs:
        sinks_mass_evolution (numpy.ndarray): Two-dimensional array containing
        the mass of each star in a star cluster at each timestep

        end_time (units.quantity): Total time of the collision

        time_step (units.quantity): Timestep of the collision

        velocity (int): Cluster velocity

        save_to (str): Path to the folder where the image should be saved

    Return:
        None
    '''

    sinks_mass_evolution = np.array(sinks_mass_evolution)
    mass_difference = sinks_mass_evolution[-1] - sinks_mass_evolution[0]

    mask = np.where(abs(mass_difference) > 1e-15)
    accreted_mass = sinks_mass_evolution - sinks_mass_evolution[0]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_facecolor('whitesmoke')

    time = np.arange(0, end_time.value_in(units.Myr),
                     time_step.value_in(units.Myr))

    plt.plot(time, accreted_mass[:, mask[0]])
    plt.grid(alpha = 0.3)
    plt.xlabel("Time [Myr]")
    plt.ylabel("Accreted mass [Msun]")
    plt.title(f"Collision with a cluster velocity = {velocity} km/s")

    if save_to is not None:
        os.makedirs(save_to, exist_ok = True)
        plt.savefig(os.path.join(save_to, f"Mass_accretion_for_cluster_velocity_{velocity}_kms.png"))
        plt.close()
    
    else:
        plt.show()

# %%

def plot_relative_mass(sinks_mass_evolution, velocity, save_to = None):
    '''
    Description:
        Plot the relative accreted mass for each star in a cluster that
        accretes mass

    Inputs:
        sinks_mass_evolution (numpy.ndarray): Two-dimensional array containing
        the mass of each star in a star cluster at each timestep

        velocity (int): Cluster velocity

        save_to (str): Path to the folder where the image should be saved

    Return:
        None
    '''

    from matplotlib.ticker import ScalarFormatter

    sinks_mass_evolution = np.array(sinks_mass_evolution)
    mass_difference = sinks_mass_evolution[-1] - sinks_mass_evolution[0]

    mask = np.where(abs(mass_difference) > 1e-15)
    relative_mass = 100*mass_difference[mask[0]]/sinks_mass_evolution[0][mask[0]]

    fig, ax = plt.subplots(figsize = (7, 5))
    ax.set_facecolor('whitesmoke')

    bins = np.logspace(np.min(np.log10(relative_mass)),
                       np.max(np.log10(relative_mass)), 30)
    plt.hist(relative_mass, bins = bins, color = '#0c2577')
    plt.axvline(5, c = "red", linestyle = "dashed", linewidth = 1.5,
                label = "Approximate limit for a second polulation")

    plt.grid(alpha = 0.3)
    plt.xscale("log")
    plt.xlabel("Relative mass difference [%]")
    ax.xaxis.set_major_formatter(ScalarFormatter())
    plt.ylabel("N")
    plt.legend()
    plt.title(f"Collision with a cluster velocity = {velocity} km/s")

    if save_to is not None:
        os.makedirs(save_to, exist_ok = True)
        plt.savefig(os.path.join(save_to, f"Relative_mass_for_cluster_velocity_{velocity}_kms.png"))
        plt.close()
    
    else:
        plt.show()

# %%

def make_3Dmap(hydro, L, N):
    '''
    Description: 
        Calculate the three-dimensional desnity of a molecular cloud

    Inputs:
        hydro (object): AMUSE hydrodynamic integrator for the molecular cloud

        L (int): Half the axis length in parsecs 

        N (int): Number of grid points

    Return: 
       rho (numpy.ndarray): Three-dimensional array for the density of the
       molecular cloud over the grid

       X (numpy.ndarray): Three-dimensional array with the X-coordinates of the
       grid's vertices

       Y (numpy.ndarray): Three-dimensional array with the Y-coordinates of the
       grid's vertices

       Z (numpy.ndarray):Three-dimensional array with the Z-coordinates of the
       grid's vertices
    '''
    # grid on the x-y-z plane
    x = np.linspace(-L, L, N + 1)
    y = np.linspace(-L, L, N + 1)
    z = np.linspace(-L, L, N + 1)
    X, Y, Z = np.meshgrid(x, y, z)

    x = X.flatten() | units.pc
    y = Y.flatten() | units.pc
    z = Z.flatten() | units.pc
    vx = 0 | units.kms
    vy = 0 | units.kms
    vz = 0 | units.kms
    # evaluate the cloud's density at the verices of the grid
    rho = hydro.get_hydro_state_at_point(x, y, z, vx, vy, vz)[0]
    rho = rho.reshape((N + 1, N + 1, N + 1))
    
    return rho, X, Y, Z




def metallicity_histogram(metallicities, N_bins,
                          collision_velocity, alpha = 1.0,
                          save_to = None):
    '''
    Description:
        This function plots a histogram of the stars' (altered) metallicity
        in log-space.
    
    Inputs:
        metallicities (np.array): An array of the metallicity values to plot
        
        N_bins (int): The number of bins
        
        alpha, optional (float): the transparency parameter
        
        collision_velocity (float): The collision velocity in kms
        (for the title)
        
        save_to (str): Path to the folder where the image should be saved
        
    Outputs:
        None
    '''
    # evenly spaced bins in log-space
    bins = np.logspace(np.min(np.log10(metallicities)),
                       np.max(np.log10(metallicities)), N_bins)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_facecolor('whitesmoke')
    ax.grid(alpha=alpha/2)

    ax.hist(metallicities, bins  = bins,
            color = '#0c2577', label =  'altered Z')
    ax.axvline(x=0.002, linestyle = '--', color = 'k', label = 'initial Z')
    ax.set_xscale('log')

    fig.supylabel('#', size = 'x-large')
    fig.supxlabel('Z', size = 'x-large')
    fig.suptitle('Collision with velocity '+ str(collision_velocity) +' kms',
                 size = 'x-large')
    ax.legend(fontsize = 14)
    plt.tight_layout()
    if save_to is not None:
        os.makedirs(save_to, exist_ok = True)
        plt.savefig(os.path.join(save_to, f"Altered_metallicities_for_cluster_velocity_{collision_velocity}_kms.png"))
        plt.close()
    
    else:
        plt.show()



def HR_diagramme(old_population, rejuvenated_population, collision_velocity,
                 x_min_inset = 3.55E+3, x_max_inset =  4E+3,
                 y_min_inset = 5E-3, y_max_inset =  5E-2,
                 alpha = 1.0, save_to = None):
    '''
    
    Description:
        This function plots an HR diagramme of two different stellar
        populations. It also contains an inseted, zoomed-in subplot.
    
    Inputs:
        old_population (object): An AMUSE particle set with attributes of 
        Luminosity and Temperature
        
        rejuvenated_population (object): An AMUSE particle set with attributes 
        of Luminosity and Temperature
        
        x_min_inset, x_max_inset optionals (floats): The range of the inset
        plot
        
        y_min_inset, y_max_inset optionals (floats): The range of the inset
        plot
        
        alpha optional (float): the transparency parameter
        
        collision_velocity (float): The collision velocity in kms
        (for the title)
        
        save_to (str): Path to the folder where the image should be saved
        
    Outputs:
        None
        
    '''
    
    fig, ax = plt.subplots( figsize=(9,6))
    ax.set_facecolor('whitesmoke')
    ax.grid(alpha=alpha/2)
    # HR diagramme of the stars with no accretion
    ax.scatter(old_population.temperature.value_in(units.K),
                    old_population.luminosity.value_in(units.LSun), 
                    c='k', s= 12,
                    label = 'Stars without accretion',zorder = 2)

    # HR diagramme of the stars with accretion
    ax.scatter(rejuvenated_population.temperature.value_in(units.K),
                      rejuvenated_population.luminosity.value_in(units.LSun), 
                      c='red', s=6,
                      label = 'Stars with accretion',zorder = 1)

    ax.set_xlim(8E+3, 3.15E+3)
    ax.set_ylim(1E-4, 5E0)
    # coloured background based on spectral type
    plt.axvspan(3150, 3500, facecolor='sandybrown', alpha=0.7, zorder = 0)
    plt.axvspan(3500, 5000, facecolor='navajowhite', alpha=0.7, zorder = 0)
    plt.axvspan(5000, 6000, facecolor='khaki', alpha=0.7, zorder = 0)
    plt.axvspan(6000, 7500, facecolor='lightyellow', alpha=0.7, zorder = 0)
    plt.axvspan(7500, 8000, facecolor='lightcyan', alpha=0.7, zorder = 0)
    #######################################
    # MS TRACK, ZOOMED IN
    ########################################

    
    # place the inset in the desired position (bounds parameter) in the graph 
    axins = ax.inset_axes(bounds =
        [.15, 0.3, .2, .4],
        xlim=(x_min_inset, x_max_inset),
        ylim=(y_min_inset, y_max_inset), xticklabels=[], yticklabels=[])

    axins.set_facecolor('whitesmoke')
    # same diagrammes, in the inset 
    axins.scatter(old_population.temperature.value_in(units.K),
                    old_population.luminosity.value_in(units.LSun), 
                    c='k', s= 12,
                    label = 'old population',zorder = 2)

    axins.scatter(rejuvenated_population.temperature.value_in(units.K),
                      rejuvenated_population.luminosity.value_in(units.LSun), 
                      c='red', s=12,
                      label = 'rejuvenated population',zorder = 1)

    # x-axis is inverted in HR diagrammes
    axins.set_xlim(x_max_inset, x_min_inset)
    axins.get_xaxis().set_visible(False)
    axins.set_ylim(y_min_inset, y_max_inset)
    axins.get_yaxis().set_visible(False)
    axins.loglog()

    ##################################################
    ax.indicate_inset_zoom(axins, edgecolor="black")
    ax.loglog()
    ax.legend(fontsize = 14)
    ax.set_xlabel("T [K]", size = 'x-large')
    ax.set_ylabel(r"L [L$_{\odot}$]", size = 'x-large')
    fig.suptitle('Collision with velocity '+ str(collision_velocity) +' kms',
                 size = 'x-large')
    plt.tight_layout()
   
    if save_to is not None:
        os.makedirs(save_to, exist_ok = True)
        plt.savefig(os.path.join(save_to, f"HR_diagramme_for_cluster_velocity_{collision_velocity}_kms.png"))
        plt.close()
    
    else:
        plt.show()


