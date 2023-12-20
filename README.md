# Simulation and Modeling in Astrophysics 2023/2024 project: YESSIR

### Authors:

	Kostas Tsalapatas - s3479765
	Erin Umuzigazuba - s3634701
	Yiqi Wu - s3805352

### General Concepts:

	This package simulates and analyses the head-on collision of a metal-poor globular cluster and a metal-rich molecular cloud. During the collision, the globular cluster accretes material from the molecular cloud.
	
	This repository can be used to verify whether this increase in mass, which is coupled with an increase in metallicity, can create a second stellar population in the cluster, as observed in most Galactic globular clusters.

### Setup & Installation:

	The repository can only be downloaded from GitHub. Download the required packages by navigating to the directory in which you downloaded the repository (make sure that `AMUSE_PROJECT_YESSIR/` is in the working directory) and run 

	```
	pip install -r requirements.txt
	```

	**Important**: The installation has only been tested on Windows (using Windows Subsystem for Linux or WSL) and Linux machines.

### Organization & Documentation:

	An overview of the repository, along with the purpose of every file and folder, can be found in `project_organization.md`. Additionally, each function includes a docstring, which can be accessed through the help() function in Python.

### Workflow:

	Our results can be reproduced as follows:
	
		- Run the `src/main_cluster_cloud_collision.py` script to simulate the collision for multiple cluster velocities between 20 and 60 km/s. For each velocity, the script generates a plot of the mass accretion evolution for each star in the cluster and a histogram of the relative accreted mass of the cluster.

		- Run the `src/updating_stellar_metallicites.py` script to update the metallicity of the stars that accreted mass, post-collision, and generate the cluster's Hertzsprung–Russell diagram. **Important**: you first need to specify which file containing the post-collision cluster must be used within the script. By default, the script uses the file containing the cluster with a velocity = 20 km/s and 1000 stars, which is the file used to make the plots in our report. 

		- Run the `src/alice_analysis.py` script to create a plot of the total mass accreted by a cluster as a function of its velocity. For a cluster velocity = 20 km/s, a plot of the mass accretion evolution and a histogram of the relative accreted mass are also created. This script uses the results obtained using the ALICE HPC cluster at Leiden University.
	
	The results of the simulations and the generated plots are saved in the `results` folder.

### Example:

	Several hours are needed to run the aforementioned scripts on a desktop computer. We therefore also provide an example script `example/cluster_cloud_collision.py` that simulates the collision for a cluster with a velocity = 20 km/s. The script also updates the metallicities of the stars after the collision and generates the cluster's Hertzsprung–Russell diagram. 

### Grade requirements:

- Minimum passing grade:
	- Main goal: Send a globular cluster, saved in a sink particle set, through a molecular cloud and observe the mass accretion. Do this by bridging the cluster's gravity code and the cloud's hydrodynamics code. Compare the results to real observations of clusters with multiple stellar populations. 
		- Done
	
	- Plot the amount of the accreted mass as a function of time for the individual stars as they passing through the cloud. (dM VS T)
		- Done

	- Plot a histogram of the final amount of accreted mass. Compare to the amount of mass accretion needed to observe multiple stellar populations. (dM VS N(dM))
		- Done

- Bonus points:
	- Repeat the experiment for different velocities and impact parameters and analyse its effect on the mass accretion.
		- Done for different velocities
	
	- Add a stellar evolution code and analyse its effect on the mass accretion. Especially the effect of stellar dynamics, e.g. wind wind-preventing accretion.
		- Done

------------

**Y**oung population

**E**mergence in old

**S**tar clusters after

**S**mashing with 

**I**nterstellar clouds leading to the

**R**egeneration of stars
