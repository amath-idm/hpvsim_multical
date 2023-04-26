'''
Calibrate HPVsim
'''

# Standard imports
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import run_sim as rs
import utils as ut

# Comment out to not run
to_run = [
    'run_calibration',
    # 'plot_calibration',
]

debug = False # Smaller runs
do_save = True

# Run settings for calibration (dependent on debug)
n_trials    = [2000, 10][debug]  # How many trials to run for calibration
n_workers   = [40, 4][debug]    # How many cores to use
storage     = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug] # Storage for calibrations


########################################################################
# Run calibration
########################################################################
def run_calib(location=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):

    sim = rs.make_sim(location)
    datafiles = ut.make_datafiles(location)[location]

    # Define the calibration parameters
    calib_pars = dict(
        beta = [0.2, 0.1, 0.3],
        # sev_dist = dict(par1=[1., 0.8, 1.2])
    )

    genotype_pars = dict(
        hpv16=dict(
            transform_prob=[8 / 1e10, 5 / 1e10, 10 / 1e10],
            sev_fn=dict(
                k = [0.2,0.15,0.25],
            ),
            dur_episomal=dict(
                par1=[0.9, .8,1],
                par2=[0.9, .8, 1],
            )
        ),
        hpv18=dict(
            transform_prob=[5 / 1e10, 3 / 1e10, 7 / 1e10],
            sev_fn=dict(
                k = [0.2,0.15,0.25],
            ),
            dur_episomal=dict(
                par1=[0.9, .8, 1],
                par2=[0.9, .8, 1],
            )
        ),
        hrhpv=dict(
            transform_prob=[4 / 1e10, 2 / 1e10, 6 / 1e10],
            sev_fn=dict(
                k = [0.2,0.15,0.25],
            ),
            dur_episomal=dict(
                par1=[0.9, .8, 1],
                par2=[0.9, .8, 1],
            )
        ),
    )

    calib = hpv.Calibration(sim, calib_pars=calib_pars, genotype_pars=genotype_pars,
                            name=f'{location}_calib',
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
def load_calib(location=None, do_plot=True, which_pars=0, save_pars=True, filestem=''):

    filename = f'{location}_calib{filestem}'
    calib = sc.load(f'results/{filename}.obj')
    if do_plot:
        sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
        sc.options(font='Libertinus Sans')
        fig = calib.plot(res_to_plot=50, plot_type='sns.boxplot', do_save=False)
        fig.suptitle(f'Calibration results, {location.capitalize()}')
        fig.tight_layout()
        fig.savefig(f'figures/{filename}_sc.png')

    if save_pars:
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=which_pars)
        sc.save(f'results/{location}_pars{filestem}_sc.obj', calib_pars)

    return calib


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    locations = ['ethiopia'] # ['drc', 'south_africa', 'kenya', 'uganda']
    filestem = '_apr26'

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        for location in locations:
            sim, calib = run_calib(location=location, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        for location in locations:
            calib = load_calib(location=location, do_plot=True, save_pars=True, filestem=filestem)
    
    T.toc('Done')