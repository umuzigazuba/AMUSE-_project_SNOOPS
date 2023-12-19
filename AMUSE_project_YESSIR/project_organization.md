# Simulation and Modeling in Astrophysics 2023/2024 project: YESSIR

Project Organization
------------

    ├── LICENSE
    │
    ├── project_organization.md
    │
    ├── __init__.py
    │
    ├── setup.py                <- Makes project pip installable (pip install -e .) so src can be imported
    │
    ├── test_environment.py     <- Checks python version
    │
    ├── requirements.txt        <- The requirements file for reproducing the analysis environment
    │
    ├── presentations           <- Slides of presentations done in class
    │   │
    │   ├── YESSIR_project_proposal.pdf           <- Presentation done before starting the project
    │   └── YESSIR_project_results.pdf            <- Presentation done when finalising the project
    │
    ├── notebooks               <- Jupyter notebooks
    │   │
    │   ├── plot_HR_proposal.ipynb                <- Create the HR diagram used for the proposal presentation
    │   └── plot_3D_collision.ipynb               <- Animate the collision in 3D
    │
    ├── src                     <- Source code for use in this project
    │   │
    │   ├── __init__.py                           <- Makes src a Python module
    │   ├── cluster_initialization.py             <- Initialize a cluster before the collision
    │   ├── molecular_cloud_initialization.py     <- Initialize a molecular cloud before the collision
    │   ├── plotters.py                           <- All plotting functions used in this project
    │   ├── utils.py                              <- Functions used during the collision
    │   ├── main_cluster_cloud_collision.py       <- Simulate the collision
    │   ├── analysis_utils.py                     <- Functions used during the analysis of our results 
    │   └── updating_stellar_metallicities.py     <- Calculate the updated metallicity of stars and create their HR diagram 
    │
    ├── example                 <- Example script
    │   │
    │   └── collision_example.py                  <- Example collision with set inital conditions/parameters
    │
    └── results                 <- Plots and data generated during our numerous runs


<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

