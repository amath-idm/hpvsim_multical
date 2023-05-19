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
import settings as set

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

    pars = dict(
        n_agents       = [10e3,1e3][debug],
        dt             = [0.25,1.0][debug],
        start          = [1960,1980][debug],
        end            = 2020,
        network        = 'default',
        genotypes      = [16, 18, 'hi5', 'ohr'],
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
        location=None, age_pyr=True, analyzers=None, debug=0, seed=0, verbose=0.2,
        do_save=False, dist_type='lognormal', marriage_scale=1, debut_bias=[0,0],
        calib_par_stem=None, calib_pars=None,
    ):

    if analyzers is None:
        analyzers = sc.autolist()
    else:
        analyzers = sc.promotetolist(analyzers)

    dflocation = location.replace(' ', '_')
    if age_pyr:
        ap = hpv.age_pyramid(
            timepoints=['2020'],
            datafile=f'data/{dflocation}_age_pyramid_reduced.csv',
            edges=np.linspace(0, 80, 9))

        analyzers += [ap]

    if calib_pars is None and calib_par_stem is not None:
        calib_pars = sc.loadobj(f'results/{location+calib_par_stem}.obj')

    # Make sim
    sim = make_sim(
        location=location,
        debug=debug,
        analyzers=analyzers,
        dist_type=dist_type,
        marriage_scale=marriage_scale,
        debut_bias=debut_bias,
        calib_pars=calib_pars
    )
    sim['rand_seed'] = seed
    sim.label = f'{location}--{seed}'

    # Run
    sim['verbose'] = verbose
    sim.run()
    sim.shrink()
        
    if do_save:
        sim.save(f'results/{dflocation}.sim')

    return sim


def run_sims(
        locations=None, age_pyr=True, debug=False, verbose=-1, analyzers=None, dist_type='lognormal',
        marriage_scale=1, debut_bias=[0,0], calib_par_stem=None, *args, **kwargs
    ):
    ''' Run multiple simulations in parallel '''
    
    kwargs = sc.mergedicts(dict(debug=debug, verbose=verbose, analyzers=analyzers, dist_type=dist_type, age_pyr=age_pyr, marriage_scale=marriage_scale, calib_par_stem=calib_par_stem, debut_bias=debut_bias), kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(location=locations), kwargs=kwargs, serial=debug, die=True)
    sims = sc.objdict({location:sim for location,sim in zip(locations, simlist)}) # Convert from a list to a dict
    
    return sims


def run_parsets(
        location=None, debug=False, verbose=-1, analyzers=None, save_results=True, **kwargs):
    ''' Run multiple simulations in parallel '''

    parsets = sc.loadobj(f'results/1_iv/{location}_pars_may18_iv_all.obj')
    kwargs = sc.mergedicts(dict(debug=debug, verbose=verbose, analyzers=analyzers), kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(calib_pars=parsets), kwargs=kwargs, serial=debug, die=True)
    msim = hpv.MultiSim(simlist)
    msim.reduce()
    if save_results:
        sc.saveobj(f'results/4_msims/{dflocation}.obj', msim.results)

    return msim


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    # locations = ['tanzania', 'uganda'] #set.locations[28:]
    # sims = run_sims(
    #     locations=locations, analyzers=[ut.dwelltime_by_genotype()],
    #     calib_par_stem='_pars_may08_sc',
    #     age_pyr=True, debug=False, verbose=.1, do_save=True)

    # location = 'mali'
    # sim = run_sim(location,  calib_par_stem='_multical_may15_pars', analyzers=[ut.dwelltime_by_genotype()], age_pyr=True, verbose=0.1, do_save=True)

    locations = ['tanzania'] #set.locations[28:]
    for location in locations:
        msim = run_parsets(location=location, save_results=True)

    T.toc('Done')

    # location = 'tanzania'
    # az = hpv.age_results(
    #     result_args=sc.objdict(
    #         hpv_prevalence=sc.objdict(
    #             years=2020,
    #             edges=[0,15,25,30,35,40,50,60,70,80,100],
    #         )
    #     )
    # )
    #
    # age_pyr = hpv.age_pyramid(
    #     timepoints=['2020'],
    #     datafile=f'data/{country}_age_pyramid_reduced.csv',
    #     edges=np.linspace(0, 100, 21))
    #
    # # calib_pars = sc.loadobj(f'results/{location}_pars_may03_sc.obj')
    # sim = make_sim(
    #     location,
    #     # calib_pars=calib_pars,
    #     analyzers=[
    #         az,
    #         age_pyr,
    #     ]
    # )
    # sim.run()
    # a = sim.get_analyzer('age_results')
    # sim.plot()
    # a.plot()
    # T.toc('Done')

