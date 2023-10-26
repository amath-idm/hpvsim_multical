'''
Utilities for multicalibration
'''

# Standard imports
import sciris as sc
import hpvsim as hpv
import hpvsim.utils as hpu
import pandas as pd
import numpy as np

import locations as set

def set_font(size=None, font='Libertinus Sans'):
    ''' Set a custom font '''
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return

def map_sb_loc(location):
    ''' Map between different representations of country names '''
    location = location.title()
    if location == "Cote Divoire": location = "Cote d'Ivoire"
    if location == "Drc": location = 'Congo Democratic Republic'
    return location


def rev_map_sb_loc(location):
    ''' Map between different representations of country names '''
    location = location.lower()
    # location = location.replace(' ', '_')
    if location == 'congo democratic republic': location = "drc"
    if location == "cote d'ivoire": location = 'cote divoire'
    return location


def make_sb_data(location=None, dist_type='lognormal', debut_bias=[0,0]):

    # Deal with missing countries and different spelling conventions
    if location in set.nosbdata_locations:
        sb_location = 'Ethiopia' # Use assumptions for Ethiopia for Sudan, South Sudan, CDI and Somalia
    else:
        sb_location = map_sb_loc(location)

    # Read in data
    sb_data_f = pd.read_csv(f'data/sb_pars_women_{dist_type}.csv')
    sb_data_m = pd.read_csv(f'data/sb_pars_men_{dist_type}.csv')

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
        f=dict(dist=distf, par1=par1f+debut_bias[0], par2=par2f),
        m=dict(dist=distm, par1=par1m+debut_bias[1], par2=par2m),
    )

    return debut


def make_datafiles(locations):
    ''' Get the relevant datafiles for the selected locations '''
    datafiles = dict()
    cancer_type_locs    = ['ethiopia', 'guinea', 'kenya', 'mali', 'mozambique', 'nigeria', 'senegal', 'south africa', 'tanzania', 'uganda', 'zimbabwe']
    cin3_type_locs      = ['guinea', 'kenya', 'nigeria', 'senegal', 'south africa', 'tanzania']
    cin1_type_locs      = ['guinea', 'kenya', 'nigeria', 'senegal', 'south africa']

    for location in locations:
        dflocation = location.replace(' ','_')
        datafiles[location] = [
            f'data/{dflocation}_cancer_cases.csv',
            f'data/{dflocation}_asr_cancer_incidence.csv',
        ]

        if location in cancer_type_locs:
            datafiles[location] += [f'data/{dflocation}_cancer_types.csv']
        # if location in cin3_type_locs:
        #     datafiles[location] += [f'data/{dflocation}_cin3_types.csv']
        # if location in cin1_type_locs:
        #     datafiles[location] += [f'data/{dflocation}_cin1_types.csv']

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

                conditions_m = ~sim.people.is_female * sim.people.alive * (sim.people.age >= (bin-1)) * (sim.people.age < bin)
                denom_inds_m = hpu.true(conditions_m)
                num_conditions_m = conditions_m * (sim.people.n_rships.sum(axis=0)>0)
                num_inds_m = hpu.true(num_conditions_m)
                self.prop_active_m[ci,bin_ind] = len(num_inds_m)/len(denom_inds_m)
        return


class dwelltime_by_genotype(hpv.Analyzer):
    '''
    Determine the age at which people with cervical cancer were causally infected and
    time spent between infection and cancer.
    '''

    def __init__(self, start_year=None, **kwargs):
        super().__init__(**kwargs)
        self.start_year = start_year
        self.years = None

    def initialize(self, sim):
        super().initialize(sim)
        self.years = sim.yearvec
        if self.start_year is None:
            self.start_year = sim['start']
        self.age_causal = dict()
        self.age_cancer = dict()
        self.dwelltime = dict()
        self.median_age_causal = dict()
        for gtype in range(sim['n_genotypes']):
            self.age_causal[gtype] = []
            self.age_cancer[gtype] = []
        for state in ['precin', 'cin', 'total']:
            self.dwelltime[state] = dict()
            for gtype in range(sim['n_genotypes']):
                self.dwelltime[state][gtype] = []

    def apply(self, sim):
        if sim.yearvec[sim.t] >= self.start_year:
            cancer_genotypes, cancer_inds = (sim.people.date_cancerous == sim.t).nonzero()
            if len(cancer_inds):
                current_age = sim.people.age[cancer_inds]
                date_exposed = sim.people.date_exposed[cancer_genotypes, cancer_inds]
                dur_precin = sim.people.dur_precin[cancer_genotypes, cancer_inds]
                dur_cin = sim.people.dur_cin[cancer_genotypes, cancer_inds]
                total_time = (sim.t - date_exposed) * sim['dt']
                for gtype in range(sim['n_genotypes']):
                    gtype_inds = hpv.true(cancer_genotypes == gtype)
                    self.dwelltime['precin'][gtype] += dur_precin[gtype_inds].tolist()
                    self.dwelltime['cin'][gtype] += dur_cin[gtype_inds].tolist()
                    self.dwelltime['total'][gtype] += total_time[gtype_inds].tolist()
                    self.age_causal[gtype] += (current_age[gtype_inds] - total_time[gtype_inds]).tolist()
                    self.age_cancer[gtype] += (current_age[gtype_inds]).tolist()
        return

    def finalize(self, sim=None):
        ''' Convert things to arrays '''
        for gtype in range(sim['n_genotypes']):
            self.median_age_causal[gtype] = np.quantile(self.age_causal[gtype], 0.5)


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
                conditions[ab] = (sim.people.age >= ab) & (sim.people.age < ab+self.binspan) & sim.people.alive & sim.people.is_female & sim.people.level0

            prop_married = sc.autolist()
            for age_cond in conditions.values():
                num_condition = age_cond & (sim.people.current_partners[0,:]>0)
                prop_married += len(hpu.true(num_condition))/len(hpv.true(age_cond))

            d = dict(age=self.bins, val=prop_married)
            df = pd.DataFrame().from_dict(d)
            df['year'] = sim.yearvec[sim.t]
            self.dfs += df

    def finalize(self, sim):
        self.df = pd.concat(self.dfs)


def lognorm_params(par1, par2):
    """
    Given the mean and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    mean = np.log(par1 ** 2 / np.sqrt(par2 ** 2 + par1 ** 2))  # Computes the mean of the underlying normal distribution
    sigma = np.sqrt(np.log(par2 ** 2 / par1 ** 2 + 1))  # Computes sigma for the underlying normal distribution

    scale = np.exp(mean)
    shape = sigma
    return shape, scale
