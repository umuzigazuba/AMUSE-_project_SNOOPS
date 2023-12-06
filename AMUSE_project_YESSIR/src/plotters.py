#from amuse.community.fi.interface import Fi

from amuse.units import units
import numpy as np
from matplotlib import pyplot as plt

# %%
def plot_snapshot_and_HR(cluster):

    '''
    Description: 
        Plot the positions of stars in a star cluster at z = 0 and its Hertzsprung-Russell diagram 

    Inputs:
        cluster (object): AMUSE particle set for the star cluster

    Return: 
        None
    '''

    v = (cluster.vx**2 + cluster.vy**2 + cluster.vz**2).sqrt()
    s = cluster.mass.value_in(units.MSun)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3))

    ax1.scatter(cluster.temperature.value_in(units.K),
                cluster.luminosity.value_in(units.LSun), 
                c = 'k',
                s = 10)
    ax1.set_xlim(8.5E+3, 2.5E+3) # Needs to depend on given cluster
    ax1.set_ylim(1.e-5, 1.e+3)
    ax1.loglog()
    ax1.set_xlabel("T [K]")
    ax1.set_ylabel("L [$L_\odot$]")

    ax2.scatter(cluster.x.value_in(units.pc), 
                cluster.y.value_in(units.pc), 
                c = v.value_in(units.kms), 
                s = s)
    plt.gca().set_aspect('equal', adjustable = 'box')
    ax2.set_xlabel("x [pc]")
    ax2.set_ylabel("y [pc]")

    plt.show()
# %%
def make_map(hydro, x_lim, y_lim, N):
    '''
    Description: 
        For a molecular cloud centered around x, y, z = 0 pc, calculate its density over a grid at z = 0 pc

    Inputs:
        hydro (object): AMUSE Fi hydrodynamic integrator for the molecular cloud

        L (int): Half the axis length in parsecs

        N (int): Number of grid points along one axis

    Return: 
        rho (numpy.ndarray): Two-dimensional array containing the density of the molecular cloud over the grid
    '''

    x = np.linspace(-x_lim, x_lim, N)
    y = np.linspace(-y_lim, y_lim, N)
    xv, yv = np.meshgrid(x, y)

    x = xv.flatten() | units.pc
    y = yv.flatten() | units.pc
    z = 0 | units.pc
    vx = 0 | units.kms
    vy = 0 | units.kms
    vz = 0 | units.kms

    rho = hydro.get_hydro_state_at_point(x, y, z, vx, vy, vz)[0]
    rho = rho.reshape((N, N))
    
    return rho

# %%  
def plot_hydro(time, hydro, x_lim, y_lim, N):
    '''
    Description: 
        Plot the log density of a molecular cloud at a given time

    Inputs: 
        time (units.quantity): Age of molecular cloud since initialization

        hydro (object): AMUSE Fi hydrodynamic integrator for the molecular cloud

        L (int): Half the axis length in parsecs

        N (int): Number of grid points along one axis

    Return:
        density_map (matplotlib.image.AxesImage): AxesImage object for the plotted log density map
    '''

    fig = plt.figure(figsize = (9, 5))
    
    rho = make_map(hydro, x_lim = x_lim, y_lim = y_lim, N = N)
    density_map = plt.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), cmap = "plasma", extent = [-x_lim, x_lim, -y_lim, y_lim])
    color_bar = fig.colorbar(density_map)
    color_bar.set_label('log density [$amu/cm^3$]', labelpad = 5)
        
    plt.title(f"Molecular cloud at time = {time.value_in(units.Myr)} Myr and z = 0 pc")
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()

    return density_map

# %%
def plot_hydro_and_star(time, hydro, star_particle, x_lim, y_lim, N, density_map_MC):
    '''
    Description: 
        Plot the log density of a molecular cloud and the position of a colliding star at a given time

    Inputs:
        time (units.quantity): Age of molecular cloud since initialization

        hydro (object): AMUSE Fi hydrodynamic integrator which contains the molecular cloud particles

        star_particle (object): AMUSE particle for the colliding star

        L (int): Half the axis length in parsecs

        N (int): Number of grid points along one axis

        density_map_MC (matplotlib.image.AxesImage): AxesImage object for the log density map of the molecular cloud before collision

    Return:
        None
    '''

    star_x = star_particle.x.value_in(units.pc)
    star_y = star_particle.y.value_in(units.pc)

    rho = make_map(hydro, x_lim = x_lim, y_lim = y_lim, N = N)
    
    fig, (ax_full, ax_zoom) = plt.subplots(nrows = 2, ncols = 1, figsize = (10.5, 7))
 
    ax_full.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), cmap = "plasma", extent = [-x_lim, x_lim, -y_lim, y_lim])
    ax_full.scatter(star_x, star_y, c = 'red')

    ax_full.set_title(f"Molecular cloud at time = {time.value_in(units.Myr)} Myr and z = 0 pc")
    ax_full.set_xlabel("x [pc]")
    ax_full.set_ylabel("y [pc]")

    ax_zoom.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), cmap = "plasma", extent = [-x_lim, x_lim, -y_lim, y_lim])
    ax_zoom.scatter(star_x, star_y, c = 'red')

    offset = x_lim/5 # parsecs

    ax_zoom.axis([star_x - offset*2, star_x + offset*2, star_y - offset, star_y + offset])
    ax_zoom.set_title(f"Zoomed in molecular cloud at time = {time.value_in(units.Myr)} Myr")
    ax_zoom.set_xlabel("x [pc]")
    ax_zoom.set_ylabel("y [pc]")
    
    colorbar_axis = fig.add_axes([0.75, 0.1, 0.02, 0.85])
    colorbar = fig.colorbar(density_map_MC, cax = colorbar_axis)
    colorbar.set_label('log density [$amu/cm^3$]', labelpad = 5)

    plt.tight_layout()
    plt.show()

# %%

def plot_cloud_and_star_cluster(time, hydro, sinks, x_lim, y_lim, N, density_map_MC):
    
    rho = make_map(hydro, x_lim = x_lim, y_lim = y_lim, N = N)

    fig = plt.figure(figsize = (9, 5))

    plt.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), cmap = "plasma", extent = [-x_lim, x_lim, -y_lim, y_lim])
    plt.scatter(sinks.position.x.value_in(units.pc), sinks.position.y.value_in(units.pc), c = "red", s = sinks.mass.value_in(units.MSun)*2)

    plt.title(f"Molecular cloud at time = {time.value_in(units.Myr)} Myr and z = 0 pc")
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.xlim([-x_lim, x_lim])
    plt.ylim([-y_lim, y_lim])

    colorbar_axis = fig.add_axes([0.95, 0.1, 0.02, 0.85])
    colorbar = plt.colorbar(density_map_MC, cax = colorbar_axis, fraction = 0.046, pad = 0.04)
    colorbar.set_label('log density [$amu/cm^3$]', labelpad = 5)

    plt.show()

# %% 
def plot_cloud_particles(time, particles_cloud):
    # Did not fill in doc string because not sure this function will stay in the final version
    '''
    Description: 

    Inputs:

    Return:
    '''

    plt.figure(figsize = (9, 5))

    plt.scatter(particles_cloud.x.value_in(units.pc), particles_cloud.y.value_in(units.pc), s = 1)
    plt.title("Molecular cloud at time = " + time.as_string_in(units.Myr))
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()

# %%

def make_3Dmap(hydro, L, N):
    '''
    Description: 
        3D density cube

    Inputs:
        hydro (object): AMUSE Fi hydrodynamic integrator

        L (int): Axis length 

        N (int): Number of grid points

    Return: 
        rho ()
    '''

    x = np.linspace(-L, L, N + 1)
    y = np.linspace(-L, L, N + 1)
    z = np.linspace(-L, L, N + 1)
    xv, yv, zv = np.meshgrid(x, y, z)

    x = xv.flatten() | units.pc
    y = yv.flatten() | units.pc
    z = zv.flatten() | units.pc
    vx = 0 | units.kms
    vy = 0 | units.kms
    vz = 0 | units.kms

    rho = hydro.get_hydro_state_at_point(x, y, z, vx, vy, vz)[0]
    rho = rho.reshape((N + 1, N + 1, N+1))
    
    return rho, xv, yv, zv
# %%


# # 3D DENSITY CUBE EXAMPLE !!LACKS A HYDRO_CODE PARTICLE SET!!

# # This is a 3d density cube
# grid = make_3Dmap(hydro_cloud, 500, 100)

# # when you take the midlle slice, the outcome is in aggreement with the plotted
# # 2D rho matrix from make_map. The middle slice of the cube correspods to
# # z = 0 (or any direction sliced for that matter)
# middle_plane = len(grid)//2

# # I am not sure if this face of the cube is indeed the x-y plane
# XY = grid[:,:,middle_plane].value_in(units.amu/units.cm**3)

# # splot the slice
# fig = plt.figure(figsize = (9, 5))
# cax = plt.imshow(np.log10(XY), extent=[-500, 500, -500, 500])

# #%%
# XZ = grid[:,middle_plane,:].value_in(units.amu/units.cm**3)
# YZ = grid[middle_plane,:,:].value_in(units.amu/units.cm**3)
