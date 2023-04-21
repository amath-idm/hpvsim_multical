'''
Calibrate HPVsim to high-burden countries and run analyses to produce estimates
of burden of cervical cancer over 2020-2060.
'''

# Additions to handle numpy multithreading
import os
os.environ.update(
    OMP_NUM_THREADS = '1',
    OPENBLAS_NUM_THREADS = '1',
    NUMEXPR_NUM_THREADS = '1',
    MKL_NUM_THREADS = '1',
)

# Standard imports
import sciris as sc
import hpvsim as hpv
import pandas as pd
import numpy as np

# Imports from this repository
import run_sim as rs
import calibration as cal
import settings as set

# Comment out to not run
to_run = [
    'run_calibration',
    # 'plot_calibration',
]


debug = False # Smaller runs
do_save = True


# Run settings for calibration (dependent on debug)
n_trials    = [3000, 2][debug]  # How many trials to run for calibration
n_workers   = [40, 1][debug]    # How many cores to use
storage     = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug] # Storage for calibrations

########################################################################
# Run calibration
########################################################################
def make_unique_priors(locations=None):
    ''' Make priors for the parameters that vary across settings '''

    unique_pars = dict()
    for location in locations:
        unique_pars[location] = dict(
            calib_pars = dict(
                beta = [0.15, 0.1, 0.25],
                # sev_dist = dict(
                #     par1 = [1.0, 0.9, 1.1]
                # ),
            ),
            genotype_pars = dict(
                hrhpv=dict(
                    transform_prob=[2 / 1e10, 1 / 1e10, 3 / 1e10],
                ),
            )
        )

    unique_pars['ethiopia']['genotype_pars']['hrhpv']['transform_prob'] = [1/1e10, 0.5/1e10, 1.5/1e10],
    unique_pars['nigeria']['genotype_pars']['hrhpv']['transform_prob'] = [2/1e10, 1.5/1e10, 3.5/1e10],
    unique_pars['senegal']['genotype_pars']['hrhpv']['transform_prob'] = [2/1e10, 1.5/1e10, 3.5/1e10],

    return unique_pars

def make_datafiles(locations):
    ''' Get the relevant datafiles for the selected locations '''
    datafiles = dict()
    cancer_type_locs    = ['ethiopia', 'guinea', 'kenya', 'mozambique', 'nigeria', 'senegal', 'south africa', 'tanzania', 'uganda']
    cin3_type_locs      = ['guinea', 'nigeria', 'senegal', 'south africa', 'tanzania']
    cin1_type_locs      = ['guinea', 'senegal', 'south africa']

    for location in locations:
        dflocation = location.replace(' ','_')
        datafiles[location] = [f'data/{dflocation}_cancer_cases.csv']

        if location in cancer_type_locs:
            datafiles[location] += [f'data/{dflocation}_cancer_types.csv']
        if location in cin3_type_locs:
            datafiles[location] += [f'data/{dflocation}_cin3_types.csv']
        if location in cin1_type_locs:
            datafiles[location] += [f'data/{dflocation}_cin1_types.csv']

    return datafiles



def run_calib(locations=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):

    # Define shared calibration parameters - same values used across sims
    common_pars = dict(
        genotype_pars=dict(
            hpv16=dict(
                transform_prob=[8/1e10, 7/1e10, 9/1e10]
            ),
            hpv18=dict(
                transform_prob=[2/1e10, 1.5/1e10, 3/1e10],
            ),
        ),
    )

    unique_pars = make_unique_priors(locations)

    sims = []
    for location in locations:
        sim = rs.make_sim(location)
        sim.label = location
        sims.append(sim)

    calib = cal.MultiCal(
        sims,
        common_pars=common_pars,
        unique_pars=unique_pars,
        name=f'multical0104',
        datafiles=make_datafiles(locations),
        load_if_exists=True,
        db_name='multical0104.db',
        total_trials=n_trials,
        n_workers=n_workers,
        storage=storage,
        keep_db=False,
    )
    calib.calibrate()

    filename = f'multical{filestem}'
    if do_plot:
        for location in locations:
            calib.plot(slabel=location, do_save=True, fig_path=f'figures/{filename}_{location}.png')
    if do_save:
        sc.saveobj(f'results/{filename}.obj', calib)

    print(f'Best pars are {calib.best_pars}')

    return sims, calib


########################################################################
# Load pre-run calibration
########################################################################
def load_calib(locations=None, do_plot=True, which_pars=0, save_pars=True):

    filename = rs.mc_filename
    calib = sc.load(f'results/{filename}.obj')
    age_bin_edges = np.array([0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 100])

    if save_pars:
        sims = []
        for location in locations:
            pars_file = f'results/{location}_{filename}_pars.obj'
            calib_pars = calib.trial_pars_to_sim_pars(slabel=location, which_pars=which_pars)
            sc.save(pars_file, calib_pars)

            # # Rerun sim with calib pars
            # az = hpv.age_results(
            #     result_args=sc.objdict(
            #         cancers=sc.objdict(
            #             years=2020,
            #             edges=age_bin_edges,
            #         )
            #     )
            # )
            # sim = rs.make_sim(location, calib_pars=calib_pars, analyzers=[az])
            # sim.run()
            # sims.append(sim)
            # sim.plot()


    if do_plot:
        sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
        sc.options(font='Libertinus Sans')
        for location in locations:
            fig = calib.plot(slabel=location, res_to_plot=20, plot_type='sns.boxplot')
            fig.suptitle(f'Calibration results, {location.capitalize()}')
            fig.tight_layout()
            fig.savefig(f'figures/{filename}_{location}.png')


    return calib, sims


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    filestem = '_apr20'

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        sims, calib = run_calib(locations=set.lo_hiv_locations, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        calib, sims = load_calib(locations=set.locations, do_plot=True, save_pars=True)

    T.toc('Done')