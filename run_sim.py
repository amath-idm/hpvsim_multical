'''
Define the HPVsim simulations for India, Nigeria, and Tanzania that are used as
the basis for the calibration, scenarios, and sweeps.

By default, all locations are run. To not run a location, comment out the line
below. For all three locations, this script should take 1-5 minutes to run.
'''

# Standard imports
import numpy as np
import sciris as sc
import hpvsim as hpv
import pandas as pd
import seaborn as sns
import pylab as pl

# Imports from this repository
import pars_data as dp
import utils as ut

#%% Settings and filepaths


mc_filename = 'multical_may01'

# Debug switch
debug = 0 # Run with smaller population sizes and in serial
do_shrink = True # Do not keep people when running sims (saves memory)

# Save settings
do_save = True
save_plots = True


#%% Simulation creation functions
def make_sim(location=None, calib_pars=None, debug=0, analyzers=[], datafile=None, seed=1,
             dist_type='lognormal', marriage_scale=1, debut_bias=[0,0]):
    ''' Define parameters, analyzers, and interventions for the simulation -- not the sim itself '''

    # # Parameters
    # location = location.replace('_', ' ')
    # if location =='cote divoire':
    #     sb_location = "cote d'ivoire"
    # else:
    #     sb_location = location

    pars = dict(
        n_agents       = [10e3,1e3][debug],
        dt             = [0.25,1.0][debug],
        start          = [1960,1980][debug],
        end            = 2020,
        network        = 'default',
        location       = location,
        debut          = ut.make_sb_data(location=location, dist_type=dist_type, debut_bias=debut_bias),
        mixing         = dp.mixing[location],
        layer_probs    = ut.make_layer_probs(location=location, marriage_scale=marriage_scale),
        partners       = dp.partners[location],
        init_hpv_dist  = dp.init_genotype_dist[location],
        init_hpv_prev  = {
            'age_brackets'  : np.array([  12,   17,   24,   34,  44,   64,    80, 150]),
            'm'             : np.array([ 0.0, 0.25, 0.6, 0.25, 0.05, 0.01, 0.0005, 0]),
            'f'             : np.array([ 0.0, 0.35, 0.7, 0.25, 0.05, 0.01, 0.0005, 0]),
        },
        ms_agent_ratio = 100,
        verbose        = 0.0,
    )

    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)

    interventions = sc.autolist()

    sim = hpv.Sim(pars=pars, interventions=interventions, analyzers=analyzers, datafile=datafile, rand_seed=seed)

    return sim



#%% Simulation running functions
def run_sim(
        location=None, analyzers=None, debug=0, seed=0, verbose=0.2,
        do_save=False, dist_type='lognormal', marriage_scale=1, debut_bias=[0,0],
    ):

    # Make sim
    sim = make_sim(
        location=location,
        debug=debug,
        analyzers=analyzers,
        dist_type=dist_type,
        marriage_scale=marriage_scale,
        debut_bias=debut_bias
    )
    sim['rand_seed'] = seed
    sim.label = f'{location}--{seed}'

    # Run
    sim['verbose'] = verbose
    sim.run()
    sim.shrink()
        
    if do_save:
        sim.save(f'results/{location}.sim')
    
    return sim


def run_sims(
        locations=None, debug=False, verbose=-1, dist_type='lognormal', marriage_scale=1, debut_bias=[0,0],
        *args, **kwargs
    ):
    ''' Run multiple simulations in parallel '''
    
    kwargs = sc.mergedicts(dict(debug=debug, verbose=verbose, dist_type=dist_type, marriage_scale=marriage_scale, debut_bias=debut_bias), kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(location=locations), kwargs=kwargs, serial=debug, die=True)
    sims = sc.objdict({location:sim for location,sim in zip(locations, simlist)}) # Convert from a list to a dict
    
    return sims


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    location = 'tanzania'
    az = hpv.age_results(
        result_args=sc.objdict(
            hpv_prevalence=sc.objdict(
                years=2020,
                edges=[0,15,25,30,35,40,50,60,70,80,100],
            )
        )
    )

    sim = make_sim(
        location,
        analyzers=[
            az
        ]
        #     ut.AFS(),
        #     ut.prop_married(),
        #     hpv.snapshot(timepoints=['2020']),
        # ],
        # dist_type='normal'
    )
    sim.run()
    a = sim.get_analyzer('age_results')
    sim.plot()
    a.plot()
    # a = sim.get_analyzer('prop_married')

    # colors = sc.gridcolors(1)
    # fig,ax = pl.subplots(figsize=(8,11))
    # sns.boxplot(data=a.df, x="age", y="val", color=colors[0], ax=ax)
    # pl.show()
    #
    #
    # # Plot age mixing
    # import matplotlib as mpl
    # import pylab as pl
    #
    # snapshot = sim.get_analyzer('snapshot')
    # people = snapshot.snapshots[0]
    # font_size = 15
    # font_family = 'Libertinus Sans'
    # pl.rcParams['font.size'] = font_size
    # pl.rcParams['font.family'] = font_family
    # fig, ax = pl.subplots(nrows=1, ncols=1, figsize=(5, 4))
    # lkey = 'm'
    # fc = people.contacts[lkey]['age_f']
    # mc = people.contacts[lkey]['age_m']
    # h = ax.hist2d(fc, mc, bins=np.linspace(0, 75, 16), density=True, norm=mpl.colors.LogNorm())
    # ax.set_xlabel('Age of female partner')
    # ax.set_ylabel('Age of male partner')
    # fig.colorbar(h[3], ax=ax)
    # ax.set_title('Marital age mixing')
    # fig.tight_layout()
    # pl.savefig(f"figures/networks.png", dpi=100)
    #
    # # Plot number of marriages, casual partnerships
    # conditions = sim.people.is_female * sim.people.alive * sim.people.level0 * (sim.people.age>=20) * (sim.people.age<25)
    # data = sim.people.current_partners[1, conditions]
    # d = np.diff(np.unique(data)).min()
    # left_of_first_bin = data.min() - float(d) / 2
    # right_of_last_bin = data.max() + float(d) / 2
    # sns.kdeplot(data)
    # pl.show()
    #
    # # Shares of women by age with 0,1,2+ current casual partners
    # binspan = 5
    # bins = np.arange(15, 50, binspan)
    # conditions = {}
    # general_conditions = sim.people.is_female * sim.people.alive * sim.people.level0
    # for ab in bins:
    #     conditions[ab] = (sim.people.age >= ab) * (sim.people.age < ab + binspan) * general_conditions
    #
    # casual_partners = {(0,1): sc.autolist(), (1,2):sc.autolist(), (2,3):sc.autolist(), (3,5):sc.autolist()}
    # for cp in casual_partners.keys():
    #     for ab,age_cond in conditions.items():
    #         this_condition = conditions[ab] * (sim.people.current_partners[1,:]>=cp[0]) * (sim.people.current_partners[1,:]<cp[1])
    #         casual_partners[cp] += len(hpv.true(this_condition))
    #
    # popsize = sc.autolist()
    # for ab, age_cond in conditions.items():
    #     popsize += len(hpv.true(age_cond))
    #
    # n_bins = len(bins)
    # partners = np.repeat([0, 1, 2, 3], n_bins)
    # allbins = np.tile(bins,4)
    # counts = np.concatenate([val for val in casual_partners.values()])
    # allpopsize = np.tile(popsize, 4)
    # shares = counts / allpopsize
    # datadict = dict(bins=allbins, partners=partners, counts=counts, popsize=allpopsize, shares=shares)
    # df = pd.DataFrame.from_dict(datadict)
    #
    # fig, ax = pl.subplots(1,1,figsize=(5, 4))
    # sns.barplot(data=df.loc[df['partners']>0], x="bins", y="shares", hue="partners",ax=ax)
    # pl.show()


    T.toc('Done')

