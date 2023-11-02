"""
Define the HPVsim simulation objects.
"""
# Additions to handle numpy multithreading
import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import pars_data as dp
import utils as ut
import locations as loc

# %% Settings and filepaths
# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Save settings
do_save = True
save_plots = True


# %% Simulation creation functions
def make_sim(location=None, calib_pars=None, debug=0, analyzers=[], datafile=None, seed=1,
             dist_type='lognormal', marriage_scale=1, debut_bias=[0, 0]):
    ''' Define parameters, analyzers, and interventions for the simulation -- not the sim itself '''

    pars = dict(
        n_agents=[10e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        start=[1960, 1980][debug],
        end=2020,
        network='default',
        genotypes=[16, 18, 'hi5', 'ohr'],
        location=location,
        debut=ut.make_sb_data(location=location, dist_type=dist_type, debut_bias=debut_bias),
        mixing=dp.mixing[location],
        layer_probs=dp.make_layer_probs(location=location, marriage_scale=marriage_scale),
        partners=dp.partners[location],
        init_hpv_dist=dp.init_genotype_dist[location],
        init_hpv_prev={
            'age_brackets': np.array([12, 17, 24, 34, 44, 64, 80, 150]),
            'm': np.array([0.0, 0.25, 0.6, 0.25, 0.05, 0.01, 0.0005, 0]),
            'f': np.array([0.0, 0.35, 0.7, 0.25, 0.05, 0.01, 0.0005, 0]),
        },
        ms_agent_ratio=100,
        verbose=0.0,
    )

    precin_pars = dict(
        hpv16=dict(dist='lognormal', par1=4, par2=20),
        hpv18=dict(dist='lognormal', par1=3, par2=15),
        hi5=dict(dist='lognormal', par1=2, par2=8),
        ohr=dict(dist='lognormal', par1=2, par2=6)
    )

    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)
        for gtype, precin_vals in precin_pars.items():
            if gtype not in pars['genotype_pars'].keys():
                pars['genotype_pars'][gtype] = dict()
            pars['genotype_pars'][gtype]['dur_precin'] = precin_vals

    interventions = sc.autolist()

    sim = hpv.Sim(pars=pars, interventions=interventions, analyzers=analyzers, datafile=datafile, rand_seed=seed)

    return sim


# %% Simulation running functions
def run_sim(
        location=None, age_pyr=True, analyzers=None, debug=0, seed=0, verbose=0.2,
        do_save=True, dist_type='lognormal', marriage_scale=1, debut_bias=[0, 0],
        calib_par_stem=None, ressubfolder=None, calib_pars=None,
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
        calib_pars = sc.loadobj(f'results/{ressubfolder}/{dflocation + calib_par_stem}.obj')

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
    # sim.shrink()

    if do_save:
        sim.save(f'results/{dflocation}.sim')

    return sim


def run_sims(
        locations=None, age_pyr=True, debug=False, verbose=-1, analyzers=None, dist_type='lognormal',
        marriage_scale=1, debut_bias=[0, 0], calib_par_stem=None, do_save=False, *args, **kwargs
):
    """ Run multiple simulations in parallel """

    kwargs = sc.mergedicts(dict(debug=debug, verbose=verbose, analyzers=analyzers, dist_type=dist_type, age_pyr=age_pyr,
                                marriage_scale=marriage_scale, calib_par_stem=calib_par_stem, debut_bias=debut_bias),
                           kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(location=locations), kwargs=kwargs, serial=debug, die=True)
    sims = sc.objdict({location: sim for location, sim in zip(locations, simlist)})  # Convert from a list to a dict

    # if do_save:
    #     for loc,sim in sims.items():
    #         sim.save(f'results/{loc}.sim')

    return sims


def run_parsets(
        location=None, debug=False, verbose=.1, analyzers=None, save_results=True, **kwargs):
    ''' Run multiple simulations in parallel '''

    dflocation = location.replace(' ', '_')
    parsets = sc.loadobj(f'results/immunovarying/{dflocation}_pars_jun15_iv_all.obj')
    kwargs = sc.mergedicts(dict(location=location, debug=debug, verbose=verbose, analyzers=analyzers), kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(calib_pars=parsets), kwargs=kwargs, serial=debug, die=True)
    msim = hpv.MultiSim(simlist)
    msim.reduce()
    if save_results:
        sc.saveobj(f'results/msims/{dflocation}.obj', msim.results)

    return msim


# %% Run as a script
if __name__ == '__main__':
    T = sc.timer()

    locations = ['uganda']  # loc.locations
    for location in locations:
        calib_pars = sc.loadobj('results/unconstrained/uganda_pars_oct26_iv.obj')
        sim = run_sim(location=location, calib_pars=calib_pars)

    # cp = sc.loadobj('results/unconstrained/tanzania_pars_oct06_iv.obj')
    # sims = run_sims(locations=locations, do_save=True)

    # alldf = sc.loadobj('results/calib_dfs_sc.obj')
    # sc_pars = dict()
    # for location in locations:
    #     thisdf = alldf[alldf.location == location]
    #     sc_pars[location] = dict(
    #         cross_imm_sev_high=thisdf.cross_imm_sev_high.iloc[0],
    #         cross_imm_sev_med=thisdf.cross_imm_sev_med.iloc[0],
    #         cross_imm_sus_high=thisdf.cross_imm_sus_high.iloc[0],
    #         cross_imm_sus_med=thisdf.cross_imm_sus_med.iloc[0],
    #         genotype_pars=dict(
    #             hi5=dict(
    #                 cin_fn=dict(k=thisdf.hi5_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50),
    #                 cancer_fn=dict(ld50=thisdf.hi5_cancer_fn_ld50.iloc[0], method='cin_integral',
    #                                **dict(k=thisdf.hi5_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50)),
    #                 dur_cin=dict(dist='lognormal', par1=thisdf.hi5_dur_cin_par1.iloc[0],
    #                              par2=thisdf.hi5_dur_cin_par2.iloc[0]),
    #                 rel_beta=thisdf.hi5_rel_beta.iloc[0]
    #             ),
    #             ohr=dict(
    #                 cin_fn=dict(k=thisdf.ohr_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50),
    #                 cancer_fn=dict(ld50=thisdf.ohr_cancer_fn_ld50.iloc[0], method='cin_integral',
    #                                **dict(k=thisdf.ohr_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50)),
    #                 dur_cin=dict(dist='lognormal', par1=thisdf.ohr_dur_cin_par1.iloc[0],
    #                              par2=thisdf.ohr_dur_cin_par2.iloc[0]),
    #                 rel_beta=thisdf.ohr_rel_beta.iloc[0]
    #             )
    #         )
    #     )

    # for location in locations:
    #
    #     sim = run_sim(location=location, calib_pars=sc_pars[location])
    #     msim = run_parsets(location=location, save_results=True)

    T.toc('Done')

