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
import utils as ut

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
                beta=[0.2, 0.1, 0.3],
                # sev_dist = dict(
                #     par1 = [1.0, 0.9, 1.1]
                # ),
            ),
            genotype_pars = dict(
                hrhpv=dict(
                    transform_prob=[3e-10, 2e-10, 5e-10],
                    sev_fn=dict(k=[0.15, 0.10, 0.2])
                ),
            )
        )

    # if 'ethiopia' in locations:
    #     unique_pars['ethiopia']['genotype_pars']['hrhpv'] = dict(
    #         transform_prob=[7 / 1e11, 5 / 1e11, 10 / 1e11],
    #         sev_fn=dict(
    #             k=[0.2, 0.15, 0.25],
    #         ),
    #     )
    #
    # if 'nigeria' in locations:
    #     unique_pars['nigeria']['genotype_pars']['hrhpv'] = dict(
    #         transform_prob=[10 / 1e11, 8 / 1e11, 12 / 1e11],
    #         sev_fn=dict(
    #             k=[0.25, 0.2, 0.3],
    #         ),
    #     )
    #
    # if 'drc' in locations:
    #     unique_pars['drc']['genotype_pars']['hrhpv'] = dict(
    #         transform_prob=[7 / 1e11, 5 / 1e11, 10 / 1e11],
    #         sev_fn=dict(
    #             k=[0.2, 0.15, 0.25],
    #         ),
    #     )
    #
    # if 'senegal' in locations:
    #     unique_pars['senegal']['genotype_pars']['hrhpv']['transform_prob'] = [2/1e10, 1.5/1e10, 3.5/1e10]

    return unique_pars



def run_calib(locations=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):

    # Define shared calibration parameters - same values used across sims
    common_pars = dict(
        genotype_pars=dict(
            hpv16=dict(transform_prob=[10e-10, 8e-10, 12e-10]),
            hpv18=dict(transform_prob=[3e-10, 2e-10, 5e-10]),
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
        datafiles=ut.make_datafiles(locations),
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
    filestem = '_may01'
    locations = ['nigeria', 'ethiopia', 'drc', 'tanzania', 'south africa', 'kenya', 'uganda']

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        sims, calib = run_calib(locations=locations, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        calib, sims = load_calib(locations=locations, do_plot=True, save_pars=True) # lo_hiv_locations

    T.toc('Done')