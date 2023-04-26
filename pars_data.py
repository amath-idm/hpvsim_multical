"""
Compilation of sexual behavior data and assumptions for three countries
with high HPV burden, for use in HPVsim.
"""


#%% Initialization

import numpy as np
import settings as set

# Initialize objects with per-country results
layer_probs          = dict()
mixing               = dict()
partners             = dict()
init_genotype_dist   = dict()


#%% LAYER PROBS
default_layer_probs = dict(
    m=np.array([
        # Share of females of each age who are married
        [0, 5,  10,    15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
        [0, 0, 0.05, 0.25, 0.70, 0.90, 0.95, 0.70, 0.75, 0.65, 0.55, 0.40, 0.40, 0.40, 0.40, 0.40], # Share of females of each age who are married
        [0, 0, 0.01, 0.01, 0.10, 0.50, 0.60, 0.70, 0.70, 0.70, 0.70, 0.80, 0.70, 0.60, 0.50, 0.60]]
    ),
    c=np.array([
        [0, 5,   10,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
        [0, 0, 0.10, 0.70, 0.80, 0.60, 0.60, 0.50, 0.20, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01], # Share of females of each age having casual relationships
        [0, 0, 0.05, 0.70, 0.80, 0.60, 0.60, 0.50, 0.50, 0.40, 0.30, 0.10, 0.05, 0.01, 0.01, 0.01]], # Share of males of each age having casual relationships
    ),
    o=np.array([
        [0, 5,   10,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
        [0, 0, 0.01, 0.05, 0.05, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01], # Share of females of each age having one-off relationships
        [0, 0, 0.01, 0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.05, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01]], # Share of males of each age having one-off relationships
    ),
)

for location in set.locations:
    layer_probs[location] = default_layer_probs

#%% INIT GENOTYPE DISTRIBUTION
default_init_genotype_dist = dict(hpv16=0.4, hpv18=0.25, hrhpv=0.35)
for location in set.locations:
    init_genotype_dist[location] = default_init_genotype_dist

# Replace for some
init_genotype_dist['ethiopia'] = dict(hpv16=0.6, hpv18=0.3, hrhpv=0.1)
init_genotype_dist['nigeria'] = dict(hpv16=0.2, hpv18=0.1, hrhpv=0.7)
init_genotype_dist['tanzania'] = dict(hpv16=0.4, hpv18=0.25, hrhpv=0.35)


#%% PARTNERS
default_partners = dict(
        m=dict(dist='poisson', par1=0.1),
        c=dict(dist='poisson', par1=0.5),
        o=dict(dist='poisson', par1=0.0),
)
for location in set.locations:
    partners[location] = default_partners


#%% MIXING
default_mixing_all = np.array([
    #       0,  5, 10, 15, 20,  25,  30,  35,  40,  45,  50,  55,  60,  65,  70,  75
    [0,     0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [5,     0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [10,    0,  0,  1,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [15,    0,  0,  1,  1,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [20,    0,  0, .5,  1,  1, .01,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [25,    0,  0,  0, .5,  1,   1, .01,   0,   0,   0,   0,   0,   0,   0,   0,   0],
    [30,    0,  0,  0,  0, .5,   1,   1, .01,   0,   0,   0,   0,   0,   0,   0,   0],
    [35,    0,  0,  0,  0, .1,  .5,   1,   1, .01,   0,   0,   0,   0,   0,   0,   0],
    [40,    0,  0,  0,  0,  0,  .1,  .5,   1,   1, .01,   0,   0,   0,   0,   0,   0],
    [45,    0,  0,  0,  0,  0,   0,  .1,  .5,   1,   1, .01,   0,   0,   0,   0,   0],
    [50,    0,  0,  0,  0,  0,   0,   0,  .1,  .5,   1,   1,  .01,  0,   0,   0,   0],
    [55,    0,  0,  0,  0,  0,   0,   0,   0,  .1,  .5,   1,   1, .01,   0,   0,   0],
    [60,    0,  0,  0,  0,  0,   0,   0,   0,   0,  .1,  .5,   1,   1, .01,   0,   0],
    [65,    0,  0,  0,  0,  0,   0,   0,   0,   0,   0,  .1,  .5,   1,   1, .01,   0],
    [70,    0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,  .1,  .5,   1,   1, .01],
    [75,    0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,  .1,  .5,   1,   1],
])

default_mixing = dict()
for k in ['m','c','o']: default_mixing[k] = default_mixing_all
for location in set.locations:
    mixing[location] = default_mixing

