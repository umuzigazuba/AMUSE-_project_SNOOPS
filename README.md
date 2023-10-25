# AMUSE_project_YESSIR

Authors: Yiqi Wu, Kostas Tsalapatas and Erin Umuzigazuba

Project: The formation of a second generation of stars in a globular cluster after a collision with a molecular cloud

General Concepts:
	We simulate the collision of a globular cluster with a certain velocity and inclination, and a dense molecular cloud that is actively forming stars, 
	in the Galactic disk. We are interested in the formation mechanism of a second generation of stars in the globular cluster. 

Experimental setup:

- Fixed parameters: mass and density for cluster and cloud

- Free parameters: impact parameters and velocity (change one at a time)

- Simulate for free fall time scale + ~2Myr years post collision
	Star formation detection
	Check if there are any new stars bound to the cluster 
	Optionally: see if mass is accreted to the first generation of stars
		Are they rejuvenated?

Codes:
- Gravity:
	BHTree
- Hydrodynamics: 
	???

Step-by-step plan:
- Initialize globular cluster and molecular cloud
- Initialize the motion of the globular cluster so that it collides with the molecular cloud 
- Simulate the collision 


Figure:
HR diagram before and after collision
Also compare first and second generations after collision
Star formation rate as a function of time

Tools to improve quality of the code/GitHub:
- ruff
- pytest
- nedbat/coveragepy
- pyscaffold => see pyscaffoldext-dsproject
- use branches! squash branches
