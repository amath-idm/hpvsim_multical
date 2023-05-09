'''
Calibrate HPVsim
'''

# Standard imports
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import run_sim as rs
import utils as ut
import settings as set

# Comment out to not run
to_run = [
    'run_calibration',
    # 'plot_calibration',
]

debug = False # Smaller runs
do_save = True

# Run settings for calibration (dependent on debug)
n_trials    = [1500, 10][debug]  # How many trials to run for calibration
n_workers   = [40, 1][debug]    # How many cores to use
storage     = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug] # Storage for calibrations


########################################################################
# Run calibration
########################################################################
def make_priors(location):
    default = dict(
            hpv16=dict(
                transform_prob=[10e-10, 8e-10, 20e-10, 1e-10],
                dur_episomal=dict(
                    par1=[2.5, 2, 5, 0.5],
                    par2=[7, 4, 10, 0.5])
            ),
            hpv18=dict(
                transform_prob=[6e-10, 4e-10, 10e-10, 1e-10],
                dur_episomal=dict(
                    par1=[2.5, 2, 3, 0.5],
                    par2=[7, 4, 10, 0.5]),
                rel_beta=[0.75, 0.7, 0.95, 0.05]
            ),
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

    genotype_pars = sc.dcp(default)

    # Old cancers: lower k and longer dur_episomal
    if location in ['tanzania']:
        genotype_pars['hpv16']=dict(
            transform_prob=[5e-10, 4e-10, 10e-10, 1e-10],
            sev_fn = dict(
                k=[0.10, 0.05, 0.2, 0.05],
            ),
            dur_episomal=dict(
                par1=[4, 3, 5, 0.5],
                par2=[12, 10, 15, 0.5])
            )
        genotype_pars['hpv18']=dict(
            sev_fn = dict(
                k=[0.10, 0.05, 0.2, 0.05],
            ),
            dur_episomal=dict(
                par1=[4, 3, 5, 0.5],
                par2=[12, 10, 15, 0.5])
        )

    return genotype_pars


def run_calib(location=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):

    sim = rs.make_sim(location)
    datafiles = ut.make_datafiles([location])[location]

    # Define the calibration parameters
    calib_pars = dict(
        beta = [0.2, 0.1, 0.3, 0.02],
        cross_imm_sus_med = [0.3, 0.2, 0.6, 0.05],
        cross_imm_sus_high = [0.5, 0.3, 0.7, 0.05],
        cross_imm_sev_med = [0.5, 0.3, 0.7, 0.05],
        cross_imm_sev_high = [0.7, 0.5, 0.9, 0.05],
        # sev_dist = dict(par1=[1.1, 1.0, 1.3])
    )
    genotype_pars = make_priors(location)

    # if location in ['mozambique', 'malawi', 'zambia', 'burundi', 'tanzania']: # Higher than 1
    #     calib_pars['sev_dist'] = dict(par1=[1.3, 1.1, 1.5, 0.005])
    # if location in ['angola', 'mali']:
    #     calib_pars['sev_dist'] = dict(par1=[1.2, 1.0, 1.4, 0.005])
    # if location in ['sudan', 'niger', 'guinea', 'zimbabwe', 'togo', 'sierra leone']:
    #     calib_pars['sev_dist'] = dict(par1=[0.9, 0.8, 1.0, 0.005])
    # if location in ['cote divoire']:
    #     calib_pars['sev_dist'] = dict(par1=[0.8, 0.7, 0.9, 0.005])

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

    fnlocation = location.replace(' ','_')
    filename = f'{fnlocation}_calib{filestem}'
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
    locations = ['nigeria', 'ethiopia', 'drc', 'south africa', 'kenya', 'angola', 'mozambique', 'ghana'] #,
    filestem = '_may08'

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        for location in locations:
            sim, calib = run_calib(location=location, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        for location in locations:
            calib = load_calib(location=location, do_plot=True, save_pars=True, filestem=filestem)
    
    T.toc('Done')