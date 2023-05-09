'''
Settings
'''

# All locations
locations = [
    'angola',       # 0
    'benin',        # 1
    'burkina faso', # 2
    'burundi',      # 3
    'cameroon',     # 4
    'chad',         # 5
    'congo',        # 6
    'cote divoire', # 7
    'drc',          # 8
    'ethiopia',     # 9
    'ghana',        # 10
    'guinea',       # 11
    'kenya',        # 12
    'madagascar',   # 13
    'malawi',       # 14
    'mali',         # 15
    'mozambique',   # 16
    'niger',        # 17
    'nigeria',      # 18
    'rwanda',       # 19
    'senegal',      # 20
    'sierra leone', # 21
    'somalia',      # 22
    'south africa', # 23
    'south sudan',  # 24
    'tanzania',     # 25
    'togo',         # 26
    'uganda',       # 27
    'zambia',       # 28
    'zimbabwe',     # 29
]

# 10 largest locations with HIV prevalence in women <2.5%
lo_hiv_locations = [
    'angola',       # 0
    'burkina faso', # 2
    'drc',          # 7
    'ethiopia',     # 8
    'madagascar',   # 12
    'mali',         # 14
    'niger',        # 16
    'nigeria',      # 17
    'senegal',      # 19
]

# Biggest 3
locations3 = [
    'nigeria',
    'ethiopia',
    'drc',
]

nosbdata_locations = ["cote d'ivoire", "cote divoire", "somalia", "south sudan"]

partitioned_locations = [
    [
        'drc',
        'south africa',
        'burkina faso',
        'senegal',
        'south sudan',
        'guinea',
        'niger',
        'ethiopia',
        'mali',
        'burundi',
    ],
    [
        'somalia',
        'zimbabwe',
        'chad',
        'cameroon',
        'kenya',
        'cote divoire',
        'congo',
        'rwanda',
        'sierra leone',
        'madagascar',
    ],
    [
        'ghana',
        'malawi',
        'nigeria',
        'zambia',
        'uganda',
        'tanzania',
        'benin',
        'angola',
        'togo',
        'mozambique',
    ]
]

fold1locations = partitioned_locations[0]+partitioned_locations[1]
fold2locations = partitioned_locations[0]+partitioned_locations[2]
fold3locations = partitioned_locations[1]+partitioned_locations[2]


code3 = [

]