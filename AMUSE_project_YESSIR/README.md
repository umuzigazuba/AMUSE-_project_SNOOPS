AMUSE project YESSIR
==============================

The formation of a second generation of stars in a globular cluster after a collision with a molecular cloud

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io

General Concepts:
------------
	We simulate the collision of a globular cluster with a certain velocity and inclination, and a dense molecular cloud that is actively forming stars, 
	in the Galactic disk. We are interested in the formation mechanism of a second generation of stars in the globular cluster. 

Experimental setup:
------------

- Fixed parameters: mass and density for cluster and cloud

- Free parameters: impact parameters and velocity (change one at a time)

- Simulate for free fall time scale + ~2Myr years post collision
	Star formation detection
	Check if there are any new stars bound to the cluster 
	Optionally: see if mass is accreted to the first generation of stars
		Are they rejuvenated?

Codes:
------------
- Gravity:
	BHTree
- Hydrodynamics: 
	???

Step-by-step plan:
------------
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


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
