'''
Utilities for multicalibration
'''

# Standard imports
import sciris as sc
import hpvsim as hpv
import pandas as pd
import numpy as np

# Comment out locations to not run
locations = [
    'angola',       # 0
    'benin',        # 1
    'burkina_faso', # 2
    'burundi',      # 3
    'cameroon',     # 4
    'chad',         # 5
    'cote_divoire', # 6
    'drc',          # 7
    'ethiopia',     # 8
    'ghana',        # 9
    'guinea',       # 10
    'kenya',        # 11
    'madagascar',   # 12
    'malawi',       # 13
    'mali',         # 14
    'mozambique',   # 15
    'niger',        # 16
    'nigeria',      # 17
    'rwanda',       # 18
    'senegal',      # 19
    'sierra_leone', # 20
    'somalia',      # 21
    'south_africa', # 22
    'south_sudan',  # 23
    'sudan',        # 24
    'tanzania',     # 25
    'togo',         # 26
    'uganda',       # 27
    'zambia',       # 28
    'zimbabwe',     # 29
]


def make_datafiles(locations):
    ''' Get the relevant datafiles for the selected locations '''
    locations = sc.promotetolist(locations)

    datafiles = dict()
    cancer_type_locs    = ['ethiopia', 'guinea', 'kenya', 'mozambique', 'nigeria', 'senegal', 'south_africa', 'tanzania', 'uganda']
    cin3_type_locs      = ['guinea', 'nigeria', 'senegal', 'south_africa', 'tanzania']
    cin1_type_locs      = ['guinea', 'senegal', 'south_africa']

    for location in locations:
        datafiles[location] = [f'data/{location}_cancer_cases.csv']

        if location in cancer_type_locs:
            datafiles[location] += [f'data/{location}_cancer_types.csv']
        if location in cin3_type_locs:
            datafiles[location] += [f'data/{location}_cin3_types.csv']
        if location in cin1_type_locs:
            datafiles[location] += [f'data/{location}_cin1_types.csv']

    return datafiles
