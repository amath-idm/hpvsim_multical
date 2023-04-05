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
storage_url     = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug] # Storage for calibrations


########################################################################
# Run calibration
########################################################################
def run_calib(location=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):

    sim = rs.make_sim(location)
    datafiles = ut.make_datafiles(location)

    # Define the calibration parameters
    calib_pars = dict(
        beta = [0.2, 0.1, 0.3],
        sev_dist = dict(par1=[1., 0.8, 1.2])
    )

    genotype_pars = dict(
        hrhpv=dict(
            transform_prob=[2 / 1e10, 1 / 1e10, 3 / 1e10],
        )
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
        calib.plot(do_save=True, fig_path=f'{ut.figfolder}/{filename}.png')
    if do_save:
        sc.saveobj(f'{ut.resfolder}/{filename}.obj', calib)

    print(f'Best pars are {calib.best_pars}')

    return sim, calib


########################################################################
# Load pre-run calibration
########################################################################
def load_calib(location=None, do_plot=True, which_pars=0, save_pars=True, do_plot_additional=False, filestem=''):

    filename = f'{location}_calib{filestem}'
    calib = sc.load(f'{ut.resfolder}/{filename}.obj')
    if do_plot:
        sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
        sc.options(font='Libertinus Sans')
        fig = calib.plot(res_to_plot=50, plot_type='sns.boxplot', do_save=True,
                         fig_path=f'{ut.figfolder}/{filename}')
        fig.suptitle(f'Calibration results, {location.capitalize()}')
        fig.tight_layout()
        fig.savefig(f'{ut.figfolder}/{filename}.png')

    if save_pars:
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=which_pars)
        sc.save(f'{ut.resfolder}/{location}_pars{filestem}.obj', calib_pars)

    if do_plot_additional:
        fig = ut.plot_trend(calib)
        ut.plot_best(calib).fig.show()

    return calib


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    locations = ['drc', 'south_africa', 'kenya', 'uganda']

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        for location in locations:
            sim, calib = run_calib(location=location, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem='_india')

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        calib = load_calib(location='nigeria', do_plot=True, save_pars=True, do_plot_additional=False, filestem='_india')

        ut.plot_calibration(locations=locations, filestem='_india')

        # calib_pars = calib.trial_pars_to_sim_pars(slabel='nigeria', which_pars=0)

        # pardict = dict(
        #     beta='beta',
        #     sev_dist_par1='sev_dist_par1',
        #     # dur_transformed_par1='dur_transformed',
        #     # hpv16_dur_episomal_par1='dur_episomal, HPV16',
        #     # hpv16_sev_fn_k='sev rate, HPV16',
        #     # hpv16_sev_fn_s='sev asymmetry, HPV16',
        #     # hpv16_sev_fn_x_infl='sev inflection, HPV16',
        #     #
        #     # hpv18_dur_episomal_par1='dur_episomal, HPV18',
        #     # hpv18_sev_fn_k='sev rate, HPV18',
        #     # hpv18_sev_fn_s='sev asymmetry, HPV18',
        #     # hpv18_sev_fn_x_infl='sev inflection, HPV18',
        #     # hpv18_transform_prob='transform prob, HPV18',
        #     #
        #     hrhpv_dur_episomal_par1='dur_episomal, hrHPV',
        #     hrhpv_sev_fn_k='sev rate, hrHPV',
        #     hrhpv_sev_fn_s='sev asymmetry, hrHPV',
        #     hrhpv_sev_fn_x_infl='sev inflection, hrHPV',
        # )

    
    T.toc('Done')