'''
Test calibration
'''

#%% Imports and settings
import sciris as sc
import hpvsim as hpv
import calibration as cal
import pylab as pl
import numpy as np

n_agents = 5e3


#%% Define the tests
def test_multical(do_plot=True):

    sc.heading('Testing multiple calibration')

    # Define individual par dicts and sims
    pars = {}
    for location in ['nigeria', 'india']:
        pars[location] = dict(n_agents=n_agents, start=1980, end=2020, dt=1., init_hpv_prev=0.8, location=location)
        age_bin_edges = np.array([0,   15,  20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70,  75,  80, 85,100])
        pars[location]['age_bin_edges'] = age_bin_edges
        pars[location]['standard_pop']  = np.array([age_bin_edges,
                                                    [0.31,.08,.08,.06,.06,.06,.06,.05,.04,.04,.03,.02,.01,.005,.005,  0,  0]])
    az = hpv.age_results(
        result_args=sc.objdict(
            cancers=sc.objdict(
                years=2020,
                edges=age_bin_edges,
            )
        )
    )

    sim_nga = hpv.Sim(pars['nigeria'], label='nigeria')
    sim_ind = hpv.Sim(pars['india'], label='india')
    sims = [sim_nga, sim_ind]

    # Define shared calibration parameters - same values used across sims
    common_pars = dict(
        calib_pars = dict(
            transm2f = [3.69, 3.5, 4],
        ),
        genotype_pars = dict(
            hpv16=dict(
                sev_fn=dict(k=[0.5, 0.2, 1.0]),
                ),
            hpv18=dict(
                sev_fn=dict(k=[0.5, 0.2, 1.0]),
            )
        )
    )

    unique_pars = dict(
        nigeria = dict(
            calib_pars = dict(
                beta = [0.05, 0.010, 0.20],
                sev_dist = dict(
                    par1=[1.1, 1, 1.4]
                ),
            ),
            genotype_pars=dict(
                hrhpv=dict(
                    sev_fn=dict(
                        k=[0.2, 0.15, 0.3],
                        x_infl=[10, 8, 20],
                    ),
                    dur_episomal=dict(
                        par1=[6, 5, 10]
                    ),
                ),
            ),
        ),
        india=dict(
            calib_pars=dict(
                beta=[0.05, 0.010, 0.20],
                sev_dist=dict(
                    par1=[1.1, 1, 1.4]
                ),
            ),
            genotype_pars=dict(
                hrhpv=dict(
                    sev_fn=dict(
                        k=[0.2, 0.15, 0.3],
                        x_infl=[10, 8, 20],
                    ),
                    dur_episomal=dict(
                        par1=[6, 5, 10]
                    ),
                ),
            ),
        ),
    )

    # Datafiles
    datafiles = dict(
        nigeria = [
            'data/nigeria_cancer_cases.csv',
            'data/nigeria_cin3_types.csv',
            'data/nigeria_cancer_types.csv',
        ],
        india = [
            'data/india_cancer_cases.csv',
            'data/india_cin1_types.csv',
            'data/india_cin3_types.csv',
            'data/india_cancer_types.csv',
        ],
    )

    calib = cal.MultiCal(
        sims,
        datafiles=datafiles,
        common_pars=common_pars,
        unique_pars=unique_pars,
        total_trials=2,
        n_workers=1,
        db_name='test_mc.db',
        keep_db=False,
    )
    calib.calibrate(die=True, tidyup=False)

    for location in ['nigeria','india']:
        calib.plot(slabel=location, res_to_plot=4, do_save=True, fig_path=f'figures/{location}_multical.png')

    # Make sure that rerunning the sims with the best pars from the calibration gives the same results
    sims = []
    for location in ['nigeria', 'india']:
        calib_pars = calib.trial_pars_to_sim_pars(slabel=location, which_pars=0)
        new_pars = sc.mergedicts(pars[location],calib_pars)
        sim = hpv.Sim(new_pars, analyzers=[az])
        sim.run()

        # Check sim results against stored results in calib
        best_run = calib.df.index[0] # Pull out the index of the best run
        year = 2020
        yind = sc.findinds(sim.results['year'], year)[0]
        calib_cancer_results = calib.age_results[location][best_run]['cancers'][2020] # Pull out the analyzer from the best run
        sim_cancer_results = sim.results['cancers_by_age'][:, yind] # Pull out the sim results from the sim run with the best pars
        az_cancer_results = sim.get_analyzer('age_results').results['cancers'][2020]

        # Do plots for visual inspection
        if do_plot:
            x = calib.age_results[location][best_run]['cancers']['bins']

            fig, ax = pl.subplots(1,1)
            ax.plot(x, sim_cancer_results, label='sim results')
            ax.plot(x, calib_cancer_results, label='calib results')
            ax.plot(x, az_cancer_results, label='analyzer results')
            ax.set_title('Cancers by age')
            ax.set_xlabel('Age')
            ax.legend()
            fig.tight_layout()
            fig.show()

        # In addition to the plots, assert that they must be equal
        assert np.allclose(calib_cancer_results,sim_cancer_results)

        sims += [sim]

    return sims, calib


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    sims, calib = test_multical(do_plot=True)

    sc.toc(T)
    print('Done.')
