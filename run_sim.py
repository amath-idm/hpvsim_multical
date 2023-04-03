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

# Locations -- comment out a line to not run
locations = [
    # 'india',        # 0
    'nigeria',      # 1
    # 'tanzania',     # 2
]
mc_filename = 'multical_apr03'

# Debug switch
debug = 0 # Run with smaller population sizes and in serial
do_shrink = True # Do not keep people when running sims (saves memory)

# Save settings
do_save = True
save_plots = True


#%% Simulation creation functions
def make_sim(location=None, calib_pars=None, debug=0, analyzers=[], datafile=None, seed=1):
    ''' Define parameters, analyzers, and interventions for the simulation -- not the sim itself '''

    # Parameters
    sbl = 'nigeria' # sexual behavior location
    pars = dict(
        n_agents       = [20e3,1e3][debug],
        dt             = [0.25,1.0][debug],
        start          = [1960,1980][debug],
        end            = 2020,
        network        = 'default',
        location       = location,
        debut          = dp.debut[sbl],
        mixing         = dp.mixing[sbl],
        layer_probs    = dp.layer_probs[sbl],
        partners       = dp.partners[sbl],
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
def run_sim(location=None, analyzers=None, calib_pars=None,
            debug=0, seed=0, verbose=0.1, end=None,
            do_save=False, die=False):

    # Make sim
    sim = make_sim(location=location, calib_pars=calib_pars, debug=debug, analyzers=analyzers, datafile=f'data/{location}_data.csv')
    sim['rand_seed'] = seed # Set seed
    sim.label = f'{location}--{seed}' # Set label

    # Run
    sim['verbose'] = verbose
    sim.run()
    sim.shrink()
        
    if do_save:
        sim.save(f'{ut.resfolder}/{location}.sim')
    
    return sim


def run_sims(locations=None, *args, **kwargs):
    ''' Run multiple simulations in parallel '''
    
    kwargs = sc.mergedicts(dict(use_calib_pars=True, debug=debug), kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(location=locations), kwargs=kwargs, serial=debug, die=True)
    sims = sc.objdict({location:sim for location,sim in zip(locations, simlist)}) # Convert from a list to a dict
    
    return sims


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    location = 'india'

    cpfiles = dict(
        default=None,
        jamie=f'results/{location}_pars_flex_sev_v6_march17.obj',
        multical=f'results/{location}_{mc_filename}_pars.obj',
    )

    which = 'multical'
    calib_pars_jc = sc.loadobj(cpfiles['jamie'])
    calib_pars_mc = sc.loadobj(cpfiles['multical'])

    az = hpv.age_results(
        result_args=sc.objdict(
            cancers_by_genotype=sc.objdict(
                years=2020,
                edges=np.array([ 0., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85.]),
            )
        )
    )

    # calib_pars_mc['genotype_pars']['hpv16']['dur_episomal'] = calib_pars_jc['genotype_pars']['hpv16']['dur_episomal']
    # calib_pars_de['genotype_pars']['hpv16']['dur_episomal'] = calib_pars_jc['genotype_pars']['hpv16']['dur_episomal']
    # calib_pars_de['genotype_pars']['hrhpv']['transform_prob'] = calib_pars_jc['genotype_pars']['hrhpv']['transform_prob']

    fig, ax = pl.subplots(1,1)
    for which,cpars in {'mc':calib_pars_mc}.items():#, 'jc':calib_pars_jc}.items():
        sim = make_sim(location=location, calib_pars=cpars, analyzers=[az])
        sim.run(verbose=0.1)
        azz = sim.get_analyzer('age_results')
        res = azz.results['cancers_by_genotype']
        age_bins = res['bins']
        for ng, gtype in enumerate(azz.glabels):
            ax.plot(age_bins, res[2020][:, ng], label=f'{gtype=}, {which=}')
        ax.legend()
    pl.show()
    # sim.plot()


    T.toc('Done')

