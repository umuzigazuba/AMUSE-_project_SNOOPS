#%%
import os
os.chdir('../src')

from time import time
from datetime import datetime
import socket
import numpy as np
import matplotlib.pyplot as plt

from molecular_cloud_initialization import *
import utils

from amuse.units import units



#%%

tot_cloud_mass = 4/3 *units.constants.pi * (15 | units.pc)**3 * ( 2.3 | units.amu * 10 / (units.cm**3))
tot_cloud_mass=tot_cloud_mass.value_in(units.MSun)
n_particles = tot_cloud_mass*15
print(tot_cloud_mass,n_particles)
star_cluster,converter_cluster = utils.make_cluster_with_vinit(20,45,207349)
#%%

init_cloud, init_cloud_converter  = make_molecular_cloud(N_cloud = n_particles,
                                                         M_cloud = tot_cloud_mass | units.MSun,
                                                         R_cloud = 15 | units.pc,
                                                         seed = 32839)

init_cloud, density_map = evolve_molecular_cloud(init_cloud, 
                                                    init_cloud_converter, 
                                                    t_end = 2 | units.Myr, 
                                                    dt = 0.2 | units.Myr, 
                                                    seed = 6829)

print("Mean density of the moelcular cloud", np.mean(init_cloud.density))

#%%
velocity = np.array([20])
t_end = np.array([3.0])
# velocity = list(reversed(velocity))
# t_end = list(reversed(t_end))

final_clusters = []

for i in range (len(velocity)):
    v = velocity[i]
    tot_t = t_end[i]
    print("Starting with cluster velocity",v)

    star_cluster,converter_cluster = utils.make_cluster_with_vinit(v,30,207349)

    stellar_evolution_code,gravhydrostellarbridge,sinks,\
        channel,particles_cloud,gravity_code,hydro_cloud = utils.hydro_gravo_stella_bridge_initialization(0.1,star_cluster,converter_cluster,\
                                                    init_cloud,init_cloud_converter)

    path = '../results/final_with_stellar/dt=0.1/'
    variable_name = str(v) + " kms"
    directory_path = os.path.join(path, str(variable_name))
    os.makedirs(directory_path, exist_ok=True)

    if __name__ == "__main__":
        start_time = time()
        date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
        print("Collision started on {}".format(socket.gethostname()),f"at time {format(date_time)}")
        
        post_collision_cluster = utils.collision_with_stellar_evolution(v,directory_path,tot_t,0.1,sinks,gravhydrostellarbridge,channel,\
                                                particles_cloud,gravity_code,hydro_cloud,stellar_evolution_code)

        print("Collision finished (running time: {0:.1f}s)".format(time() - start_time))
#%%
print(star_cluster)
#%%
#retrieve the data and plot trends 

names = [20,30,40,50,60]
mass_diff = []
for i in range(len(names)):
    name = names[i]
    file_path = os.path.join('../results/final_without_stellar//dt=0.2/all txt/' ,f"{name}.txt")
    data = np.loadtxt(file_path)
    x_values = np.arange(1, len(data) + 1)
    
    diff = data[-1] - data[0]
    mask = np.where(diff > 1e-15)
    dm = np.sum(diff)
    mass_diff.append(dm)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.scatter(names,mass_diff)
ax.plot(names,mass_diff)
ax.set_xticks(names)

plt.xlabel("Cluster velocity [km/s]")
plt.ylabel("Total accreted mass [MSun]")
plt.title("Total accreted mass VS Velocity")
plt.savefig('../results/final_without_stellar//dt=0.2/mass_velocity.png')


#%%
plt.plot(x_values, data)
plt.xlabel('Time Period')
plt.ylabel('Accretion')
plt.title('Accretion Trend')
plt.grid(True)
plt.show()
plt.close()

sinks_mass_snapshots = data
relative_mass = data - data[0,:]
mass_difference = sinks_mass_snapshots[-1] - sinks_mass_snapshots[0]

# Remove the column with the most increase
plt.plot(x_values, relative_mass)
plt.xlabel('Time (0.1 Myr)')
plt.ylabel('Accretion')
plt.title("Collision with velocity 60 kms")
plt.grid(True)
plt.savefig('../results/final_without_stellar/60 kms/accretion_over_time_60kms.png')
plt.show()
plt.close()

mask = np.where(mass_difference> 1e-15)


mass_ratio = np.array(mass_difference)[mask[0]]/np.array(sinks_mass_snapshots[0])[mask[0]]
plt.hist(mass_ratio*100, bins  = 30)
plt.xlabel("Relative mass difference [%]")
plt.title('Collision with velocity 60 kms')
plt.savefig('../results/final_without_stellar/60 kms/accretion_hist_60kms.png')
plt.show()
plt.close()

# %%
gravity_code.stop()
hydro_cloud.stop()
#gravhydrobridge.stop()

# %%
