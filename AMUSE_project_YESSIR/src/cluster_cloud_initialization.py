#%%
from amuse.units import units
from amuse.lab import Particles, new_kroupa_mass_distribution
from amuse.community.seba.interface import SeBa
import numpy as np
from amuse.ic.kingmodel import new_king_model
import matplotlib.pyplot as plt
from amuse.ext.masc import new_star_cluster
# %%

def make_globular_cluster():
    