'''
Utilities for multicalibration
'''

# Standard imports
import sciris as sc
import hpvsim as hpv
import hpvsim.utils as hpu
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

nosbdata_locations = ["cote d'ivoire", "somalia", "south sudan", "sudan"]

def map_sb_loc(location):
    ''' Map between different representations of country names '''
    # location = location.replace('_', ' ')
    # if location =='cote divoire': location = "cote d'ivoire"
    location = location.title()
    if location == "Drc": location = 'Congo Democratic Republic'
    return location


def rev_map_sb_loc(location):
    ''' Map between different representations of country names '''
    location = location.lower()
    location = location.replace(' ', '_')
    if location == 'congo_democratic_republic': location = "drc"
    if location == "cote_d'ivoire": location = 'cote_divoire'
    return location


def make_sb_data(location=None):

    # Deal with missing countries and different spelling conventions
    if location in nosbdata_locations:
        sb_location = 'Ethiopia' # Use assumptions for Ethiopia for Sudan, South Sudan, CDI and Somalia
    else:
        sb_location = map_sb_loc(location)

    # Read in data
    sb_data_f = pd.read_csv('data/sb_pars_women.csv')
    sb_data_m = pd.read_csv('data/sb_pars_men.csv')

    try:
        distf = sb_data_f.loc[sb_data_f["location"]==sb_location,"dist"].iloc[0]
        par1f = sb_data_f.loc[sb_data_f["location"]==sb_location,"par1"].iloc[0]
        par2f = sb_data_f.loc[sb_data_f["location"]==sb_location,"par2"].iloc[0]
        distm = sb_data_m.loc[sb_data_m["location"]==sb_location,"dist"].iloc[0]
        par1m = sb_data_m.loc[sb_data_m["location"]==sb_location,"par1"].iloc[0]
        par2m = sb_data_m.loc[sb_data_m["location"]==sb_location,"par2"].iloc[0]
    except:
        print(f'No data for {sb_location=}, {location=}')

    debut = dict(
        f=dict(dist=distf, par1=par1f, par2=par2f),
        m=dict(dist=distm, par1=par1m, par2=par2m),
    )

    return debut



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


class AFS(hpv.Analyzer):
    def __init__(self, bins=None, cohort_starts=None, **kwargs):
        super().__init__(**kwargs)
        self.bins = bins or np.arange(12,31,1)
        self.cohort_starts = cohort_starts
        self.binspan = self.bins[-1]-self.bins[0]

    def initialize(self, sim):
        super().initialize()
        if self.cohort_starts is None:
            first_cohort = sim['start'] + sim['burnin'] - 5
            last_cohort = sim['end']-self.binspan
            self.cohort_starts = sc.inclusiverange(first_cohort, last_cohort)
            self.cohort_ends = self.cohort_starts+self.binspan
            self.n_cohorts = len(self.cohort_starts)
            self.cohort_years = np.array([sc.inclusiverange(i,i+self.binspan) for i in self.cohort_starts])

        self.prop_active_f = np.zeros((self.n_cohorts,self.binspan+1))
        self.prop_active_m = np.zeros((self.n_cohorts,self.binspan+1))

    def apply(self, sim):
        if sim.yearvec[sim.t] in self.cohort_years:
            cohort_inds, bin_inds = sc.findinds(self.cohort_years, sim.yearvec[sim.t])
            for ci,cohort_ind in enumerate(cohort_inds):
                bin_ind = bin_inds[ci]
                bin = self.bins[bin_ind]

                conditions_f = sim.people.is_female * sim.people.alive * (sim.people.age >= (bin-1)) * (sim.people.age < bin) * sim.people.level0
                denom_inds_f = hpu.true(conditions_f)
                num_conditions_f = conditions_f * (sim.people.n_rships.sum(axis=0)>0)
                num_inds_f = hpu.true(num_conditions_f)
                self.prop_active_f[cohort_ind,bin_ind] = len(num_inds_f)/len(denom_inds_f)

                conditions_m = ~sim.people.is_female * sim.people.alive * (sim.people.age >= (bin-1)) * (sim.people.age < bin) * sim.people.level0
                denom_inds_m = hpu.true(conditions_m)
                num_conditions_m = conditions_m * (sim.people.n_rships.sum(axis=0)>0)
                num_inds_m = hpu.true(num_conditions_m)
                self.prop_active_m[ci,bin_ind] = len(num_inds_m)/len(denom_inds_m)
        return


class prop_married(hpv.Analyzer):
    def __init__(self, bins=None, years=None, includelast=True, yearstride=5, binspan=5, **kwargs):
        super().__init__(**kwargs)
        self.bins = bins or np.arange(15,50,binspan)
        self.years = years
        self.dfs = sc.autolist()
        self.df = None
        self.includelast = includelast
        self.yearstride = yearstride
        self.binspan = binspan

    def initialize(self, sim):
        super().initialize()
        if self.years is None:
            start = sim['start'] + sim['burnin']
            end = sim['end']
            self.years = np.arange(start, end, self.yearstride)
            if self.includelast:
                if end not in self.years:
                    self.years = np.append(self.years, end)

    def apply(self, sim):
        if sim.yearvec[sim.t] in self.years:

            conditions = dict()
            for ab in self.bins:
                conditions[ab] = (sim.people.age >= ab) & (sim.people.age < ab+self.binspan)

            prop_married = sc.autolist()
            for age_cond in conditions.values():
                num_condition = age_cond & (sim.people.current_partners[0,:]>0)
                prop_married += len(hpv.true(num_condition))/len(hpv.true(age_cond))

            d = dict(age=self.bins, val=prop_married)
            df = pd.DataFrame().from_dict(d)
            df['year'] = sim.yearvec[sim.t]
            self.dfs += df

    def finalize(self, sim):
        self.df = pd.concat(self.dfs)
