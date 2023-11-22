#from amuse.community.fi.interface import Fi

from amuse.units import units
import numpy as np
from matplotlib import pyplot as plt

#%%

# plot molecular cloud density function (smooth)
def make_map(hydro, L, N):
    '''
    Description: 2D density map at the z = 0 plane (HardCoded)

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
    '''
    Description: ploting a 2D density map using make_map function
    '''
    fig = plt.figure(figsize = (9, 5))
    
    rho = make_map(hydro, L = L, N = N)
    cax = plt.imshow(np.log10(rho.value_in(units.amu/units.cm**3)), extent=[-L, L, -L, L]) # , vmin = 0, vmax = 5
    cbar = fig.colorbar(cax)
    cbar.set_label('log density [$amu/cm^3$]', labelpad = 5)
        
    plt.title("Molecular cloud at time = " + time.as_string_in(units.Myr))
    plt.xlabel("x [pc]")
    plt.ylabel("x [pc]")
    plt.show()


#%%

def make_3Dmap(hydro, L, N):
    '''
    Description: 3D density cube

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
    
    return rho
#%%


# 3D DENSITY CUBE EXAMPLE !!LACKS A HYDRO_CODE PARTICLE SET!!

# This is a 3d density cube
grid = make_3Dmap(hydro_cloud, 500, 100)

# when you take the midlle slice, the outcome is in aggreement with the plotted
# 2D rho matrix from make_map. The middle slice of the cube correspods to
# z = 0 (or any direction sliced for that matter)
middle_plane = len(grid)//2

# I am not sure if this face of the cube is indeed the x-y plane
XY = grid[:,:,middle_plane].value_in(units.amu/units.cm**3)

# splot the slice
fig = plt.figure(figsize = (9, 5))
cax = plt.imshow(np.log10(XY), extent=[-500, 500, -500, 500])

#%%
XZ = grid[:,middle_plane,:].value_in(units.amu/units.cm**3)
YZ = grid[middle_plane,:,:].value_in(units.amu/units.cm**3)
