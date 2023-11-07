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
cal_type = ['unconstrained', 'immunovarying'][0]  # Whether to run the unconstrained or immunovarying calibration
debug = False  # If True, this will do smaller runs that can be run locally for debugging
do_save = True

# Run settings for calibration (dependent on debug)
n_trials = [3000, 10][debug]  # How many trials to run for calibration
n_workers = [40, 1][debug]  # How many cores to use
storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations


########################################################################
# Run calibration
########################################################################
def make_priors(add_1618=True):
    default = dict(
        cin_fn=dict(k=[.25, .2, .4, 0.01]),
        rel_beta=[0.9, 0.8, 1.2, 0.05]
    )

    genotype_pars = dict(
        hi5=sc.dcp(default),
        ohr=sc.dcp(default)
    )

    if add_1618:
        hpv16 = dict(
            cin_fn=dict(k=[.3, .2, .4, 0.01]),
        )
        hpv18 = dict(
            cin_fn=dict(k=[.3, .2, .4, 0.01]),
            rel_beta=[0.9, 0.8, 1.2, 0.05]
        )
        genotype_pars['hpv16'] = hpv16
        genotype_pars['hpv18'] = hpv18

    return genotype_pars


def run_calib(location=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem='', mc_gpars=None):

    sim = rs.make_sim(location, calib_pars=mc_gpars)
    datafiles = ut.make_datafiles([location])[location]

    # Define the calibration parameters
    calib_pars = dict(
        beta=[0.2, 0.1, 0.34, 0.02],
        m_cross_layer=[0.3, 0.1, 0.7, 0.05],
        m_partners=dict(
            c=dict(par1=[0.2, 0.1, 0.6, 0.02])
        ),
        f_cross_layer=[0.1, 0.05, 0.5, 0.05],
        f_partners=dict(
            c=dict(par1=[0.2, 0.1, 0.6, 0.02])
        )
    )

    if mc_gpars is None: add_1618 = True
    else: add_1618 = False
    genotype_pars = make_priors(add_1618=add_1618)

    # Save some extra sim results
    extra_sim_result_keys = ['cancers', 'cancer_incidence', 'asr_cancer_incidence']

    calib = hpv.Calibration(sim, calib_pars=calib_pars, genotype_pars=genotype_pars,
                            name=f'{location}_calib_final',
                            datafiles=datafiles,
                            extra_sim_result_keys=extra_sim_result_keys,
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
    locations = ['cote divoire']  # loc.locations[:10]
    filestem = '_nov06'

    if cal_type == 'immunovarying':
        mc_gpars = dict(
            genotype_pars=dict(
                hpv16=dict(
                    cin_fn=dict(k=[0.3, 0.1, 0.5, 0.01]),
                ),
                hpv18=dict(
                    cin_fn=dict(k=[0.3, 0.1, 0.5, 0.01]),
                    rel_beta=0.95
                ),
            )
        )

    elif cal_type == 'unconstrained':
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
