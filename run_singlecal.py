"""
This file is used to run both the immuno-varying and the unconstrained calibrations as described in
the HPVsim multicalibration paper (reference forthcoming).

Instructions: Go to the CONFIGURATIONS section on lines 29-36 to set up the script before running it.
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
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import run_sim as rs
import utils as ut
import locations as loc

# CONFIGURATIONS TO BE SET BY USERS BEFORE RUNNING
to_run = [
    'run_calibration',  # Make sure this is uncommented if you want to _run_ the calibrations (usually on VMs)
    # 'plot_calibration',  # Make sure this is uncommented if you want to _plot_ the calibrations (usually locally)
]
cal_type = ['unconstrained', 'immunovarying'][1]  # Whether to run the unconstrained or immunovarying calibration
debug = False  # If True, this will do smaller runs that can be run locally for debugging
do_save = True

# Run settings for calibration (dependent on debug)
n_trials = [1200, 10][debug]  # How many trials to run for calibration
n_workers = [40, 1][debug]  # How many cores to use
storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations


########################################################################
# Run calibration
########################################################################
def make_priors(add_1618=True):
    default = dict(
        hi5=dict(
            transform_prob=[3e-10, 2e-10, 5e-10, 1e-10],
            sev_fn=dict(k=[0.05, 0.04, 0.8, 0.01]),
            dur_episomal=dict(
                par1=[2.5, 2, 3, 0.5],
                par2=[7, 4, 10, 0.5]),
            rel_beta=[0.75, 0.7, 1.25, 0.05]
        ),
        ohr=dict(
            transform_prob=[3e-10, 2e-10, 5e-10, 1e-10],
            sev_fn=dict(k=[0.05, 0.04, 0.8, 0.01]),
            dur_episomal=dict(
                par1=[2.5, 2, 3, 0.5],
                par2=[7, 4, 10, 0.5]),
            rel_beta=[0.75, 0.7, 1.25, 0.05]
        ),
    )

    if add_1618:
        default['hpv16'] = dict(
            transform_prob=[10e-10, 4e-10, 20e-10, 1e-10],
            sev_fn=dict(
                k=[0.25, 0.15, 0.4, 0.05],
            ),
            dur_episomal=dict(
                par1=[2.5, 1.5, 5, 0.5],
                par2=[7, 4, 15, 0.5])
        )
        default['hpv18'] = dict(
            transform_prob=[6e-10, 4e-10, 10e-10, 1e-10],
            sev_fn=dict(
                k=[0.2, 0.1, 0.35, 0.05],
            ),
            dur_episomal=dict(
                par1=[2.5, 1.5, 3, 0.5],
                par2=[7, 4, 15, 0.5]),
            rel_beta=[0.75, 0.7, 0.95, 0.05]
        )

    genotype_pars = sc.dcp(default)

    return genotype_pars


def run_calib(location=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem='', mc_gpars=None):

    sim = rs.make_sim(location, calib_pars=mc_gpars)
    datafiles = ut.make_datafiles([location])[location]

    # Define the calibration parameters
    calib_pars = dict(
        beta=[0.2, 0.1, 0.3, 0.02],
        cross_imm_sus_med=[0.3, 0.2, 0.6, 0.05],
        cross_imm_sus_high=[0.5, 0.3, 0.7, 0.05],
        cross_imm_sev_med=[0.5, 0.3, 0.7, 0.05],
        cross_imm_sev_high=[0.7, 0.5, 0.9, 0.05],
        sev_dist=dict(par1=[2.0, 1.5, 5.0, 0.05])
    )
    if location in ['south africa', 'kenya']:
        calib_pars['beta'] = [0.2, 0.14, 0.3, 0.02]
        calib_pars['cross_imm_sev_high'] = [0.5, 0.3, 0.5, 0.05]
        calib_pars['cross_imm_sus_high'] = [0.5, 0.3, 0.7, 0.05]
        calib_pars['sev_dist'] = dict(par1=[3.0, 1.0, 5.0, 0.1])
    if location == 'tanzania':
        calib_pars['beta'] = [0.2, 0.14, 0.3, 0.02]
        calib_pars['cross_imm_sus_med'] = [0.25, 0.2, 0.3, 0.05]
        calib_pars['cross_imm_sus_high'] = [0.4, 0.3, 0.5, 0.05]
        calib_pars['cross_imm_sev_med'] = [0.4, 0.3, 0.5, 0.05]
        calib_pars['cross_imm_sev_high'] = [0.5, 0.4, 0.6, 0.05]
        calib_pars['sev_dist'] = dict(par1=[1.0, 0.6, 2.0, 0.1])

    if mc_gpars is None: add_1618 = True
    else: add_1618 = False
    genotype_pars = make_priors(add_1618=add_1618)

    calib = hpv.Calibration(sim, calib_pars=calib_pars, genotype_pars=genotype_pars,
                            name=f'{location}_calib_final',
                            datafiles=datafiles,
                            total_trials=n_trials, n_workers=n_workers,
                            storage=storage
                            )
    calib.calibrate()
    filename = f'{location}_calib{filestem}'
    if do_plot:
        calib.plot(do_save=True, fig_path=f'figures/{filename}.png')
    if do_save:
        sc.saveobj(f'results/{filename}.obj', calib)

    print(f'Best pars are {calib.best_pars}')

    return sim, calib


########################################################################
# Load pre-run calibration
########################################################################
def load_calib(location=None, do_plot=True, which_pars=0, save_pars=True, filestem='', ressubfolder=None,
               figsubfolder=None):
    fnlocation = location.replace(' ', '_')
    filename = f'{fnlocation}_calib{filestem}'
    calib = sc.load(f'results/{ressubfolder}/{filename}.obj')
    if do_plot:
        sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
        sc.options(font='Libertinus Sans')
        fig = calib.plot(res_to_plot=200, plot_type='sns.boxplot', do_save=False)
        fig.suptitle(f'Calibration results, {location.capitalize()}')
        fig.tight_layout()
        fig.savefig(f'figures/{figsubfolder}/{filename}_sc.png')

    if save_pars:
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=which_pars)
        trial_pars = sc.autolist()
        for i in range(100):
            trial_pars += calib.trial_pars_to_sim_pars(which_pars=i)
        sc.save(f'results/{ressubfolder}/{location}_pars{filestem}_iv.obj', calib_pars)
        sc.save(f'results/{ressubfolder}/{location}_pars{filestem}_iv_all.obj', trial_pars)

    return calib


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()
#    rererun_locations = ['uganda', 'zambia']  # [1.0, 1.15, 1.2, 1.1, 1.05, 1.2, 1.1]
    rerererun_locations = ['tanzania']  # [1.0, 1.15, 1.2, 1.1, 1.05, 1.2, 1.1]
    locations = rerererun_locations
    filestem = '_jun15'
    # ressubfolder = '1a_iv'
    # figsubfolder = '6_may19iv'

    if cal_type == 'immunovarying':
        mc_gpars = dict(
            genotype_pars=dict(
                hpv16=dict(
                    transform_prob=1.3e-9,
                    sev_fn=dict(form='logf2', k=0.15, x_infl=0, ttc=30),
                    dur_episomal=dict(dist='lognormal', par1=1.5, par2=4)
                ),
                hpv18=dict(
                    transform_prob=9e-10,
                    sev_fn=dict(form='logf2', k=0.1, x_infl=0, ttc=30),
                    dur_episomal=dict(dist='lognormal', par1=2, par2=12),
                    rel_beta=0.95
                ),
            )
        )

    elif cal_type == 'unconstrained' in to_run:
        mc_gpars = None
    else:
        raise ValueError('Need to define which calibration to run.')

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        for location in locations:
            sim, calib = run_calib(mc_gpars=mc_gpars, location=location, n_trials=n_trials, n_workers=n_workers,
                                   do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        for location in locations:
            calib = load_calib(location=location, do_plot=True, save_pars=True, filestem=filestem,
                               ressubfolder=cal_type, figsubfolder=cal_type)

    T.toc('Done')
