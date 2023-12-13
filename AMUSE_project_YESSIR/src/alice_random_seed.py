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

rng = np.random.default_rng(12345)



tot_cloud_mass = 4/3 *units.constants.pi * ( 15 | units.pc)**3 * ( 2.3 | units.amu * 10 / (units.cm**3))
tot_cloud_mass=tot_cloud_mass.value_in(units.MSun)
n_particles = tot_cloud_mass*15
print(tot_cloud_mass,n_particles)

velocity = np.array([20,30,40,50,60])
t_end = np.array([2.8,2.0,1.6,1.2,1.0])
velocity = list(reversed(velocity))
t_end = list(reversed(t_end))


for i in range(20):
    cloud_seed = []
    cluster_seed = []
    random1 = rng.integers(low = 1, high = 999999, size=1)
    random2= rng.integers(low = 1, high = 999999, size=1)

    cloud_seed.append(random1)
    cluster_seed.append(random2)

    path = '../results/alice/random_seeds/'

    init_cloud, init_cloud_converter  = make_molecular_cloud(N_cloud = n_particles,
                                                            M_cloud = tot_cloud_mass | units.MSun,
                                                            R_cloud = 15 | units.pc,
                                                            seed = random1)

    init_cloud, density_map = evolve_molecular_cloud(init_cloud, 
                                                        init_cloud_converter, 
                                                        t_end = 2 | units.Myr, 
                                                        dt = 0.2 | units.Myr, 
                                                        seed = random1)

    print("Mean density of the moelcular cloud", np.mean(init_cloud.density))


    for k in range (len(velocity)):
        v = velocity[i]
        tot_t = t_end[i]
        print("Starting with cluster velocity",v)

        star_cluster,converter_cluster = utils.make_cluster_with_vinit(v,30,random2,200)

        gravhydrobridge,sinks,channel,particles_cloud,gravity_code,\
            hydro_cloud = utils.AMUSE_bridge_initialization(0.1,star_cluster,converter_cluster,\
                                                        init_cloud,init_cloud_converter)
        
        variable_name = "seed "+str(i)
        directory_path = os.path.join(path, str(variable_name))
        os.makedirs(directory_path, exist_ok=True)

        if __name__ == "__main__":
            start_time = time()
            date_time = datetime.utcfromtimestamp(start_time).strftime('%Y-%m-%d_%H-%M-%S')
            print("Collision started on {}".format(socket.gethostname()),f"at time {format(date_time)}")
            
            post_collision_cluster = utils.let_them_collide_and_save(v,directory_path,tot_t,0.1,sinks,gravhydrobridge,channel,\
                                                    particles_cloud,gravity_code,hydro_cloud)

            print("Collision finished (running time: {0:.1f}s)".format(time() - start_time))


