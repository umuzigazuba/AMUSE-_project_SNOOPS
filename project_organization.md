# Simulation and Modeling in Astrophysics 2023/2024 project: YESSIR

Project Organization
------------

    ├── LICENSE
    │
    ├── README.md
    │
    ├── project_organization.md
    │
    ├── __init__.py             <- Makes AMUSE_PROJECT_YESSIR a Python module
    │
    ├── setup.py                <- Makes project pip installable (pip install -e .) so src can be imported
    │
    ├── test_environment.py     <- Checks python version
    │
    ├── requirements.txt        <- List of the packages that need to be installed in order to use our package
    │
    ├── AMUSE_YESSIR_report.pdf <- Report
    │
    ├── presentations           <- Slides of the presentations done in class
    │   │
    │   ├── YESSIR_project_proposal.pdf           <- Presentation done before starting the project
    │   └── YESSIR_project_results.pdf            <- Presentation done when finalising the project
    │
    ├── notebooks               <- Jupyter notebooks
    │   │
    │   ├── plot_HR_proposal.ipynb                <- Create the HR diagram shown in the proposal presentation
    │   └── plot_3D_collision.ipynb               <- Animate the collision in 3D
    │
    ├── src                     <- Source code used in this project
    │   │
    │   ├── __init__.py                           <- Makes src a Python module
    │   ├── cluster_initialization.py             <- Initialize a cluster before the collision
    │   ├── molecular_cloud_initialization.py     <- Initialize a molecular cloud before the collision
    │   ├── plotters.py                           <- Plotting functions used in this project
    │   ├── utils.py                              <- Functions used during the collision
    │   ├── main_cluster_cloud_collision.py       <- Simulate the collision
    │   ├── analysis_utils.py                     <- Functions used during the analysis of our results 
    │   ├── updating_stellar_metallicities.py     <- Calculate the updated metallicity of stars and create their HR diagram 
    │   └── alice_anlysis.py                      <- Analysis of the results exclusively obtained using ALICE
    │
    ├── example                 <- Example script
    │   │
    │   ├── collision_example.py                  <- Example collision with set inital conditions/parameters
    │   └── results_example                       <- Results obtained using the collision_example.py script
    │
    └── results                 <- Plots and data generated during our numerous runs
        │
        ├── alice                                 <- Results exclusively obtained using ALICE: convergence test, simulations of 
        │                                            clusters with 200, 1000 and 10000 stars and simulations of clusters with
        │                                            200 stars and ten different velocities, repeated using 19 different random seeds
        └── final_with_stellar_evolution          <- Results obtained using our personal computers: simulations of clusters with
                                                     1000 stars and five different velocities. Includes snapshots of the molecular cloud at
                                                     z = 0 pc and the cluster at each timestep


<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

