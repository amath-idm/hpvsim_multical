"""
Compilation of sexual behavior data and assumptions for 30 SSA countries
"""


#%% Initialization

import numpy as np
import locations as loc
import utils as ut
import pandas as pd
import sciris as sc

# Initialize objects with per-country results
layer_probs = dict()
mixing = dict()
partners = dict()
init_genotype_dist = dict()


# %% LAYER PROBS
default_layer_probs = dict(
    m=np.array([
        # Share of people of each age who are married
        [0, 5,  10,    15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
        [0, 0, 0.05, 0.25, 0.70, 0.90, 0.95, 0.70, 0.75, 0.65, 0.55, 0.40, 0.40, 0.40, 0.40, 0.40],  # Females
        [0, 0, 0.01, 0.01, 0.10, 0.50, 0.60, 0.70, 0.70, 0.70, 0.70, 0.80, 0.70, 0.60, 0.50, 0.60]]  # Males
    ),
    c=np.array([
        # Share of people of each age in casual partnerships
        [0, 5,   10,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
        [0, 0, 0.10, 0.70, 0.80, 0.60, 0.60, 0.50, 0.20, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],  # Females
        [0, 0, 0.05, 0.70, 0.80, 0.60, 0.60, 0.50, 0.50, 0.40, 0.30, 0.10, 0.05, 0.01, 0.01, 0.01]],  # Males
    ),
)

def make_layer_probs(location=None, marriage_scale=1):
    # Deal with missing countries and different spelling conventions
    if location in loc.nosbdata_locations:
        sb_location = 'Ethiopia' # Use assumptions for Ethiopia for Sudan, South Sudan, CDI and Somalia
    else:
        sb_location = ut.map_sb_loc(location)

    # Read in data and write to layer_probs
    prop_married = pd.read_csv(f'data/prop_married.csv')
    vals = np.array(prop_married.loc[prop_married["Country"] == sb_location, ["15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"]])[0]
    layer_probs = sc.dcp(default_layer_probs)
    layer_probs['m'][1][3:10] = vals/100
    layer_probs['m'][1]*=marriage_scale

    # Make individual country adjustments
    if location == 'angola':
        layer_probs['m'][1] *= .7
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.15, 0.05, 0.01, 0.01, 0.01
        ])
    if location == 'benin':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.25, 0.1, 0.05, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'burkina faso':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.4, 0.5, 0.4, 0.3, 0.2, 0.2, 0.1, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01
        ])

    # if location=='benin':
    #     layer_probs['m'][1][9:] *= .15
    if location=='burundi':
        layer_probs['m'][1][4:]*=.7
    if location=='cameroon':
        layer_probs['m'][1]*=.7
    if location=='chad':
        layer_probs['m'][1]*=.7
        layer_probs['m'][1][9:]*=.15
    if location=='congo':
        layer_probs['m'][1]*=.5
    if location=='cote divoire':
        layer_probs['m'][1]*=.7
    if location=='drc':
        layer_probs['m'][1]*=.7
        # layer_probs['m'][1][9:] *= .15
    if location=='ethiopia':
        layer_probs['m'][1]*=.7
        # layer_probs['m'][1][9:]*=.15
    if location=='ghana':
        layer_probs['m'][1]*=.6
    if location in ['guinea', 'nigeria', 'senegal', 'sierra leone', 'sudan']:
        layer_probs['m'][1]*=.7
    if location in ['zambia']:
        layer_probs['m'][1]*=.9
    if location in ['uganda', 'zimbabwe']:
        layer_probs['m'][1]*=.7
    if location in ['kenya', 'madagascar', 'malawi', 'mozambique', 'rwanda', 'togo']:
        layer_probs['m'][1]*=.5
        # layer_probs['m'][1][9:]*=.15
    if location=='tanzania':
        layer_probs['m'][1]*=.7
        # layer_probs['m'][1][9:]=.35
        # layer_probs['m'][1][9:]*=.15
    if location=='south africa':
        layer_probs['m'][1]*=.8
        # layer_probs['m'][1][7:]*=.5


    # Increasing numbers of casual partners for women in Tanzania as they get older
    if location == 'burundi':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.3, 0.2, 0.2, 0.2, 0.3, 0.25, 0.2, 0.1, 0.05, 0.01, 0.01, 0.01
        ])
    if location == 'cameroon':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.5, 0.6, 0.4, 0.4, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'chad':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.2, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'congo':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.3, 0.4, 0.5, 0.7, 0.6, 0.6, 0.4, 0.2, 0.1, 0.05, 0.05, 0.05
        ])
    if location == 'cote divoire':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.5, 0.3, 0.1, 0.05, 0.05, 0.05, 0.05
        ])
    if location == 'drc':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.5, 0.4, 0.3, 0.3, 0.2, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'ethiopia':
        layer_probs['c'][1] = 3*np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.07, 0.1, 0.1, 0.1, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'ghana':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.65, 0.6, 0.5, 0.5, 0.5, 0.3, 0.3, 0.1, 0.05, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'guinea':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.3, 0.3, 0.4, 0.5, 0.6, 0.6, 0.5, 0.3, 0.1, 0.01, 0.01, 0.01
        ])
    if location == 'kenya':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.3, 0.3, 0.4, 0.5, 0.6, 0.6, 0.6, 0.6, 0.1, 0.01, 0.01, 0.01
        ])
    if location == 'madagascar':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.4, 0.4, 0.5, 0.5, 0.5, 0.3, 0.2, 0.1, 0.1, 0.05, 0.01, 0.01, 0.01
        ])
    if location == 'malawi':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.4, 0.5, 0.4, 0.3, 0.2, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'mali':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.45, 0.2, 0.1, 0.1, 0.1, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'mozambique':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.7, 0.7, 0.6, 0.4, 0.4, 0.25, 0.25, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'nigeria':
        layer_probs['c'][1] = 2*np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.3, 0.4, 0.3, 0.2, 0.2, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'rwanda':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.6, 0.3, 0.2, 0.2, 0.2, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'senegal':
        layer_probs['c'][1] = 7*np.array([
            # 0, 5,   10,   15,   20,  25,  30,    35,   40,   45,   50,  55,  60,  65,   70,   75
            0,   0, 0.04, 0.05, 0.04, 0.05, 0.05, 0.1, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'sierra leone':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.8, 0.7, 0.6, 0.6, 0.5, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'somalia':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.5, 0.6, 0.7, 0.8, 0.6, 0.2, 0.1, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'south africa':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.8, 0.8, 0.65, 0.7, 0.6, 0.3, 0.3, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01
        ])
    # if location == 'benin':
    #     layer_probs['c'][1] = np.array([
    #         # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
    #         0,   0, 0.1, 0.7, 0.8, 0.6, 0.6, 0.5, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
    #     ])
    if location == 'tanzania':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,    25,  30,  35,  40,   45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.25, 0.3, 0.3, 0.5, 0.5, 0.5, 0.6, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2
            # 0,   0, 0.1, 0.25, 0.5, 0.75, 0.6, 0.6, 0.7, 0.75, 0.3, 0.3, 0.1, 0.1, 0.05, 0.01
        ])
        layer_probs['c'][2] = layer_probs['c'][1]
    if location == 'togo':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.65, 0.5, 0.2, 0.2, 0.2, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'uganda':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.7, 0.5, 0.3, 0.3, 0.2, 0.2, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'zambia':
        layer_probs['c'][1] = np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.6, 0.4, 0.4, 0.3, 0.5, 0.6, 0.5, 0.4, 0.1, 0.01, 0.01, 0.01, 0.01
        ])
    if location == 'zimbabwe':
        layer_probs['c'][1] = 2*np.array([
            # 0, 5,  10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75
            0,   0, 0.1, 0.15, 0.2, 0.3, 0.3, 0.25, 0.25, 0.2, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01
        ])

    return layer_probs





#%% INIT GENOTYPE DISTRIBUTION
default_init_genotype_dist = dict(hpv16=0.4, hpv18=0.25, hi5=0.25, ohr=.1)
for location in loc.locations:
    init_genotype_dist[location] = default_init_genotype_dist


#%% PARTNERS
default_partners = dict(
        m=dict(dist='poisson', par1=0.1),
        c=dict(dist='poisson', par1=0.5),
)
for location in loc.locations:
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
for k in ['m', 'c']: default_mixing[k] = default_mixing_all
for location in loc.locations:
    mixing[location] = default_mixing
