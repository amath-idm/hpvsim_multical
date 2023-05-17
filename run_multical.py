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
    # 'run_calibration',
    'plot_calibration',
]


debug = False # Smaller runs
do_save = True


# Run settings for calibration (dependent on debug)
n_trials    = [10000, 2][debug]  # How many trials to run for calibration
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
                beta=[0.2, 0.1, 0.3, 0.01],
                cross_imm_sus_med=[0.3, 0.2, 0.6, 0.05],
                cross_imm_sus_high=[0.5, 0.3, 0.7, 0.05],
                cross_imm_sev_med=[0.5, 0.3, 0.7, 0.05],
                cross_imm_sev_high=[0.7, 0.5, 0.9, 0.05],
            ),
            genotype_pars = dict(
                hi5=dict(
                    transform_prob=[4e-10, 2e-10, 6e-10, 1e-10],
                    sev_fn=dict(k=[0.15, 0.05, 0.2, 0.01]),
                    dur_episomal=dict(
                        par1=[2.5, 2, 3, 0.5],
                        par2=[7, 4, 10, 0.5]),
                    rel_beta=[0.75, 0.7, 1.25, 0.05]
                ),
                ohr=dict(
                    transform_prob=[4e-10, 2e-10, 6e-10, 1e-10],
                    sev_fn=dict(k=[0.15, 0.05, 0.2, 0.01]),
                    dur_episomal=dict(
                        par1=[2.5, 2, 3, 0.5],
                        par2=[7, 4, 10, 0.5]),
                    rel_beta=[0.75, 0.7, 1.25, 0.05]
                ),
            )
        )

    return unique_pars



def run_calib(locations=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):

    # Define shared calibration parameters - same values used across sims
    common_pars = dict(
        genotype_pars=dict(
            hpv16=dict(
                transform_prob=[10e-10, 4e-10, 20e-10, 1e-10],
                sev_fn=dict(
                    k=[0.25, 0.1, 0.4, 0.05],
                ),
                dur_episomal=dict(
                    par1=[2.5, 2, 5, 0.5],
                    par2=[7, 4, 15, 0.5])
            ),
            hpv18=dict(
                transform_prob=[6e-10, 4e-10, 10e-10, 1e-10],
                sev_fn=dict(
                    k=[0.25, 0.1, 0.4, 0.05],
                ),
                dur_episomal=dict(
                    par1=[2.5, 2, 3, 0.5],
                    par2=[7, 4, 15, 0.5]),
                rel_beta=[0.75, 0.7, 0.95, 0.05]
            ),
        )
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
def load_calib(filestem=None, locations=None, do_plot=True, which_pars=0, save_pars=True):

    calib = sc.load(f'results/multical{filestem}.obj')

    if save_pars:
        sims = []
        for location in locations:
            pars_file = f'results/{location}_multical{filestem}_pars.obj'
            calib_pars = calib.trial_pars_to_sim_pars(slabel=location, which_pars=which_pars)
            sc.save(pars_file, calib_pars)

    if do_plot:
        sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
        sc.options(font='Libertinus Sans')
        for location in locations:
            fig = calib.plot(slabel=location, res_to_plot=50, plot_type='sns.boxplot')
            fig.suptitle(f'Calibration results, {location.capitalize()}')
            fig.tight_layout()
            fig.savefig(f'figures/7_may15mc/multical{filestem}_{location}.png')


    return calib, sims


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    filestem = '_may15'
    locations = set.locations #['nigeria', 'ethiopia', 'drc', 'tanzania', 'south africa', 'kenya', 'uganda', 'angola', 'mozambique', 'ghana']

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        sims, calib = run_calib(locations=locations, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        calib, sims = load_calib(filestem=filestem, locations=locations, do_plot=True, save_pars=True) # lo_hiv_locations

    T.toc('Done')