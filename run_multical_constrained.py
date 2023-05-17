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
n_trials    = [4000, 2][debug]  # How many trials to run for calibration
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
                sev_dist = dict(par1=[1.0, 0.5, 1.5]),
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
              do_plot=False, do_save=True, filestem='', mc_gpars=None):

    unique_pars = make_unique_priors(locations)

    sims = []
    for location in locations:
        sim = rs.make_sim(location, calib_pars=mc_gpars)
        sim.label = location
        sims.append(sim)

    calib = cal.MultiCal(
        sims,
        # common_pars=common_pars,
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
            fig = calib.plot(slabel=location, res_to_plot=20, plot_type='sns.boxplot')
            fig.suptitle(f'Calibration results, {location.capitalize()}')
            fig.tight_layout()
            fig.savefig(f'figures/multical{filestem}_{location}.png')


    return calib, sims


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    filestem = '_may17'
    locations = set.locations #['nigeria', 'ethiopia', 'drc', 'tanzania', 'south africa', 'kenya', 'uganda', 'angola', 'mozambique', 'ghana']
    mc_gpars = dict(
        genotype_pars=dict(
            hpv16=dict(
                transform_prob=5e-10,
                sev_fn=dict(k=0.1),
                dur_episomal=dict(par1=2, par2=11)
            ),
            hpv18=dict(
                transform_prob=7e-10,
                sev_fn=dict(k=0.15),
                dur_episomal=dict(par1=2, par2=15),
                rel_beta=0.85
            ),
        )
    )

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        sims, calib = run_calib(mc_gpars=mc_gpars, locations=locations, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        calib, sims = load_calib(filestem=filestem, locations=locations, do_plot=True, save_pars=True) # lo_hiv_locations

    T.toc('Done')

    # best_pars = {'hpv16_transform_prob': 5e-10, 'hpv16_sev_fn_k': 0.1,
    #                                       'hpv16_dur_episomal_par1': 2.0, 'hpv16_dur_episomal_par2': 11.0,
    #                                       'hpv18_transform_prob': 7.000000000000001e-10,
    #                                       'hpv18_sev_fn_k': 0.15000000000000002, 'hpv18_dur_episomal_par1': 2.0,
    #                                       'hpv18_dur_episomal_par2': 15.0, 'hpv18_rel_beta': 0.85, 'angola_beta': 0.19,
    #                                       'angola_cross_imm_sus_med': 0.45, 'angola_cross_imm_sus_high': 0.7,
    #                                       'angola_cross_imm_sev_med': 0.5, 'angola_cross_imm_sev_high': 0.75,
    #                                       'angola_hi5_transform_prob': 3e-10, 'angola_hi5_sev_fn_k': 0.09,
    #                                       'angola_hi5_dur_episomal_par1': 2.0, 'angola_hi5_dur_episomal_par2': 4.0,
    #                                       'angola_hi5_rel_beta': 1.0, 'angola_ohr_transform_prob': 3e-10,
    #                                       'angola_ohr_sev_fn_k': 0.15000000000000002,
    #                                       'angola_ohr_dur_episomal_par1': 2.5, 'angola_ohr_dur_episomal_par2': 10.0,
    #                                       'angola_ohr_rel_beta': 0.7, 'benin_beta': 0.13,
    #                                       'benin_cross_imm_sus_med': 0.55, 'benin_cross_imm_sus_high': 0.3,
    #                                       'benin_cross_imm_sev_med': 0.65, 'benin_cross_imm_sev_high': 0.7,
    #                                       'benin_hi5_transform_prob': 6e-10, 'benin_hi5_sev_fn_k': 0.14,
    #                                       'benin_hi5_dur_episomal_par1': 2.5, 'benin_hi5_dur_episomal_par2': 8.5,
    #                                       'benin_hi5_rel_beta': 0.8999999999999999, 'benin_ohr_transform_prob': 3e-10,
    #                                       'benin_ohr_sev_fn_k': 0.1, 'benin_ohr_dur_episomal_par1': 2.0,
    #                                       'benin_ohr_dur_episomal_par2': 6.5, 'benin_ohr_rel_beta': 0.75,
    #                                       'burkina faso_beta': 0.12000000000000001,
    #                                       'burkina faso_cross_imm_sus_med': 0.5, 'burkina faso_cross_imm_sus_high': 0.5,
    #                                       'burkina faso_cross_imm_sev_med': 0.4, 'burkina faso_cross_imm_sev_high': 0.7,
    #                                       'burkina faso_hi5_transform_prob': 2e-10, 'burkina faso_hi5_sev_fn_k': 0.08,
    #                                       'burkina faso_hi5_dur_episomal_par1': 2.5,
    #                                       'burkina faso_hi5_dur_episomal_par2': 9.0, 'burkina faso_hi5_rel_beta': 1.25,
    #                                       'burkina faso_ohr_transform_prob': 5e-10, 'burkina faso_ohr_sev_fn_k': 0.13,
    #                                       'burkina faso_ohr_dur_episomal_par1': 2.5,
    #                                       'burkina faso_ohr_dur_episomal_par2': 8.0, 'burkina faso_ohr_rel_beta': 0.75,
    #                                       'burundi_beta': 0.24000000000000002, 'burundi_cross_imm_sus_med': 0.45,
    #                                       'burundi_cross_imm_sus_high': 0.4, 'burundi_cross_imm_sev_med': 0.3,
    #                                       'burundi_cross_imm_sev_high': 0.9, 'burundi_hi5_transform_prob': 2e-10,
    #                                       'burundi_hi5_sev_fn_k': 0.12000000000000001,
    #                                       'burundi_hi5_dur_episomal_par1': 2.0, 'burundi_hi5_dur_episomal_par2': 4.0,
    #                                       'burundi_hi5_rel_beta': 1.2, 'burundi_ohr_transform_prob': 4e-10,
    #                                       'burundi_ohr_sev_fn_k': 0.060000000000000005,
    #                                       'burundi_ohr_dur_episomal_par1': 2.5, 'burundi_ohr_dur_episomal_par2': 6.5,
    #                                       'burundi_ohr_rel_beta': 1.05, 'cameroon_beta': 0.12000000000000001,
    #                                       'cameroon_cross_imm_sus_med': 0.5, 'cameroon_cross_imm_sus_high': 0.65,
    #                                       'cameroon_cross_imm_sev_med': 0.65, 'cameroon_cross_imm_sev_high': 0.75,
    #                                       'cameroon_hi5_transform_prob': 4e-10,
    #                                       'cameroon_hi5_sev_fn_k': 0.12000000000000001,
    #                                       'cameroon_hi5_dur_episomal_par1': 2.0, 'cameroon_hi5_dur_episomal_par2': 9.5,
    #                                       'cameroon_hi5_rel_beta': 0.8999999999999999,
    #                                       'cameroon_ohr_transform_prob': 4e-10, 'cameroon_ohr_sev_fn_k': 0.09,
    #                                       'cameroon_ohr_dur_episomal_par1': 3.0, 'cameroon_ohr_dur_episomal_par2': 9.0,
    #                                       'cameroon_ohr_rel_beta': 1.25, 'chad_beta': 0.18,
    #                                       'chad_cross_imm_sus_med': 0.6, 'chad_cross_imm_sus_high': 0.3,
    #                                       'chad_cross_imm_sev_med': 0.7, 'chad_cross_imm_sev_high': 0.65,
    #                                       'chad_hi5_transform_prob': 4e-10, 'chad_hi5_sev_fn_k': 0.09,
    #                                       'chad_hi5_dur_episomal_par1': 2.5, 'chad_hi5_dur_episomal_par2': 9.0,
    #                                       'chad_hi5_rel_beta': 0.95, 'chad_ohr_transform_prob': 5e-10,
    #                                       'chad_ohr_sev_fn_k': 0.19, 'chad_ohr_dur_episomal_par1': 2.5,
    #                                       'chad_ohr_dur_episomal_par2': 5.0, 'chad_ohr_rel_beta': 1.0,
    #                                       'congo_beta': 0.22, 'congo_cross_imm_sus_med': 0.5,
    #                                       'congo_cross_imm_sus_high': 0.3, 'congo_cross_imm_sev_med': 0.35,
    #                                       'congo_cross_imm_sev_high': 0.75, 'congo_hi5_transform_prob': 2e-10,
    #                                       'congo_hi5_sev_fn_k': 0.08, 'congo_hi5_dur_episomal_par1': 2.0,
    #                                       'congo_hi5_dur_episomal_par2': 6.0, 'congo_hi5_rel_beta': 0.95,
    #                                       'congo_ohr_transform_prob': 4e-10, 'congo_ohr_sev_fn_k': 0.07,
    #                                       'congo_ohr_dur_episomal_par1': 2.5, 'congo_ohr_dur_episomal_par2': 5.5,
    #                                       'congo_ohr_rel_beta': 1.1, 'cote divoire_beta': 0.1,
    #                                       'cote divoire_cross_imm_sus_med': 0.4, 'cote divoire_cross_imm_sus_high': 0.5,
    #                                       'cote divoire_cross_imm_sev_med': 0.55,
    #                                       'cote divoire_cross_imm_sev_high': 0.9,
    #                                       'cote divoire_hi5_transform_prob': 5e-10,
    #                                       'cote divoire_hi5_sev_fn_k': 0.060000000000000005,
    #                                       'cote divoire_hi5_dur_episomal_par1': 2.0,
    #                                       'cote divoire_hi5_dur_episomal_par2': 7.0, 'cote divoire_hi5_rel_beta': 1.25,
    #                                       'cote divoire_ohr_transform_prob': 6e-10, 'cote divoire_ohr_sev_fn_k': 0.19,
    #                                       'cote divoire_ohr_dur_episomal_par1': 2.0,
    #                                       'cote divoire_ohr_dur_episomal_par2': 5.5, 'cote divoire_ohr_rel_beta': 0.7,
    #                                       'drc_beta': 0.1, 'drc_cross_imm_sus_med': 0.30000000000000004,
    #                                       'drc_cross_imm_sus_high': 0.7, 'drc_cross_imm_sev_med': 0.55,
    #                                       'drc_cross_imm_sev_high': 0.6, 'drc_hi5_transform_prob': 2e-10,
    #                                       'drc_hi5_sev_fn_k': 0.09, 'drc_hi5_dur_episomal_par1': 2.5,
    #                                       'drc_hi5_dur_episomal_par2': 9.5, 'drc_hi5_rel_beta': 0.8999999999999999,
    #                                       'drc_ohr_transform_prob': 2e-10, 'drc_ohr_sev_fn_k': 0.08,
    #                                       'drc_ohr_dur_episomal_par1': 2.5, 'drc_ohr_dur_episomal_par2': 10.0,
    #                                       'drc_ohr_rel_beta': 1.05, 'ethiopia_beta': 0.15000000000000002,
    #                                       'ethiopia_cross_imm_sus_med': 0.6, 'ethiopia_cross_imm_sus_high': 0.55,
    #                                       'ethiopia_cross_imm_sev_med': 0.7, 'ethiopia_cross_imm_sev_high': 0.75,
    #                                       'ethiopia_hi5_transform_prob': 3e-10, 'ethiopia_hi5_sev_fn_k': 0.2,
    #                                       'ethiopia_hi5_dur_episomal_par1': 2.0, 'ethiopia_hi5_dur_episomal_par2': 8.5,
    #                                       'ethiopia_hi5_rel_beta': 1.15, 'ethiopia_ohr_transform_prob': 6e-10,
    #                                       'ethiopia_ohr_sev_fn_k': 0.08, 'ethiopia_ohr_dur_episomal_par1': 2.0,
    #                                       'ethiopia_ohr_dur_episomal_par2': 4.0, 'ethiopia_ohr_rel_beta': 0.85,
    #                                       'ghana_beta': 0.28, 'ghana_cross_imm_sus_med': 0.4,
    #                                       'ghana_cross_imm_sus_high': 0.6000000000000001,
    #                                       'ghana_cross_imm_sev_med': 0.35, 'ghana_cross_imm_sev_high': 0.55,
    #                                       'ghana_hi5_transform_prob': 2e-10, 'ghana_hi5_sev_fn_k': 0.12000000000000001,
    #                                       'ghana_hi5_dur_episomal_par1': 2.5, 'ghana_hi5_dur_episomal_par2': 9.0,
    #                                       'ghana_hi5_rel_beta': 1.1, 'ghana_ohr_transform_prob': 2e-10,
    #                                       'ghana_ohr_sev_fn_k': 0.08, 'ghana_ohr_dur_episomal_par1': 2.5,
    #                                       'ghana_ohr_dur_episomal_par2': 10.0, 'ghana_ohr_rel_beta': 0.95,
    #                                       'guinea_beta': 0.11, 'guinea_cross_imm_sus_med': 0.25,
    #                                       'guinea_cross_imm_sus_high': 0.45, 'guinea_cross_imm_sev_med': 0.5,
    #                                       'guinea_cross_imm_sev_high': 0.55, 'guinea_hi5_transform_prob': 6e-10,
    #                                       'guinea_hi5_sev_fn_k': 0.05, 'guinea_hi5_dur_episomal_par1': 3.0,
    #                                       'guinea_hi5_dur_episomal_par2': 10.0, 'guinea_hi5_rel_beta': 0.75,
    #                                       'guinea_ohr_transform_prob': 5e-10, 'guinea_ohr_sev_fn_k': 0.19,
    #                                       'guinea_ohr_dur_episomal_par1': 2.0, 'guinea_ohr_dur_episomal_par2': 8.5,
    #                                       'guinea_ohr_rel_beta': 1.2, 'kenya_beta': 0.11,
    #                                       'kenya_cross_imm_sus_med': 0.55, 'kenya_cross_imm_sus_high': 0.5,
    #                                       'kenya_cross_imm_sev_med': 0.6000000000000001,
    #                                       'kenya_cross_imm_sev_high': 0.8, 'kenya_hi5_transform_prob': 5e-10,
    #                                       'kenya_hi5_sev_fn_k': 0.05, 'kenya_hi5_dur_episomal_par1': 2.5,
    #                                       'kenya_hi5_dur_episomal_par2': 9.5, 'kenya_hi5_rel_beta': 1.25,
    #                                       'kenya_ohr_transform_prob': 4e-10, 'kenya_ohr_sev_fn_k': 0.15000000000000002,
    #                                       'kenya_ohr_dur_episomal_par1': 2.5, 'kenya_ohr_dur_episomal_par2': 4.5,
    #                                       'kenya_ohr_rel_beta': 0.7, 'madagascar_beta': 0.11,
    #                                       'madagascar_cross_imm_sus_med': 0.30000000000000004,
    #                                       'madagascar_cross_imm_sus_high': 0.45, 'madagascar_cross_imm_sev_med': 0.65,
    #                                       'madagascar_cross_imm_sev_high': 0.7, 'madagascar_hi5_transform_prob': 4e-10,
    #                                       'madagascar_hi5_sev_fn_k': 0.08, 'madagascar_hi5_dur_episomal_par1': 3.0,
    #                                       'madagascar_hi5_dur_episomal_par2': 9.0, 'madagascar_hi5_rel_beta': 0.75,
    #                                       'madagascar_ohr_transform_prob': 6e-10, 'madagascar_ohr_sev_fn_k': 0.07,
    #                                       'madagascar_ohr_dur_episomal_par1': 2.0,
    #                                       'madagascar_ohr_dur_episomal_par2': 9.0, 'madagascar_ohr_rel_beta': 0.85,
    #                                       'malawi_beta': 0.11, 'malawi_cross_imm_sus_med': 0.30000000000000004,
    #                                       'malawi_cross_imm_sus_high': 0.65,
    #                                       'malawi_cross_imm_sev_med': 0.6000000000000001,
    #                                       'malawi_cross_imm_sev_high': 0.6, 'malawi_hi5_transform_prob': 6e-10,
    #                                       'malawi_hi5_sev_fn_k': 0.07, 'malawi_hi5_dur_episomal_par1': 2.5,
    #                                       'malawi_hi5_dur_episomal_par2': 5.0, 'malawi_hi5_rel_beta': 0.7,
    #                                       'malawi_ohr_transform_prob': 4e-10, 'malawi_ohr_sev_fn_k': 0.07,
    #                                       'malawi_ohr_dur_episomal_par1': 2.0, 'malawi_ohr_dur_episomal_par2': 5.5,
    #                                       'malawi_ohr_rel_beta': 1.2, 'mali_beta': 0.25, 'mali_cross_imm_sus_med': 0.4,
    #                                       'mali_cross_imm_sus_high': 0.65, 'mali_cross_imm_sev_med': 0.35,
    #                                       'mali_cross_imm_sev_high': 0.65, 'mali_hi5_transform_prob': 5e-10,
    #                                       'mali_hi5_sev_fn_k': 0.14, 'mali_hi5_dur_episomal_par1': 2.5,
    #                                       'mali_hi5_dur_episomal_par2': 9.0, 'mali_hi5_rel_beta': 0.7999999999999999,
    #                                       'mali_ohr_transform_prob': 3e-10, 'mali_ohr_sev_fn_k': 0.05,
    #                                       'mali_ohr_dur_episomal_par1': 2.5, 'mali_ohr_dur_episomal_par2': 8.5,
    #                                       'mali_ohr_rel_beta': 0.7, 'mozambique_beta': 0.11,
    #                                       'mozambique_cross_imm_sus_med': 0.4, 'mozambique_cross_imm_sus_high': 0.65,
    #                                       'mozambique_cross_imm_sev_med': 0.55, 'mozambique_cross_imm_sev_high': 0.55,
    #                                       'mozambique_hi5_transform_prob': 3e-10, 'mozambique_hi5_sev_fn_k': 0.07,
    #                                       'mozambique_hi5_dur_episomal_par1': 2.5,
    #                                       'mozambique_hi5_dur_episomal_par2': 5.5, 'mozambique_hi5_rel_beta': 0.85,
    #                                       'mozambique_ohr_transform_prob': 6e-10, 'mozambique_ohr_sev_fn_k': 0.13,
    #                                       'mozambique_ohr_dur_episomal_par1': 2.5,
    #                                       'mozambique_ohr_dur_episomal_par2': 8.0, 'mozambique_ohr_rel_beta': 0.75,
    #                                       'niger_beta': 0.15000000000000002,
    #                                       'niger_cross_imm_sus_med': 0.30000000000000004,
    #                                       'niger_cross_imm_sus_high': 0.35, 'niger_cross_imm_sev_med': 0.4,
    #                                       'niger_cross_imm_sev_high': 0.65, 'niger_hi5_transform_prob': 6e-10,
    #                                       'niger_hi5_sev_fn_k': 0.13, 'niger_hi5_dur_episomal_par1': 2.0,
    #                                       'niger_hi5_dur_episomal_par2': 9.0, 'niger_hi5_rel_beta': 1.0,
    #                                       'niger_ohr_transform_prob': 6e-10, 'niger_ohr_sev_fn_k': 0.08,
    #                                       'niger_ohr_dur_episomal_par1': 2.5, 'niger_ohr_dur_episomal_par2': 10.0,
    #                                       'niger_ohr_rel_beta': 0.85, 'nigeria_beta': 0.11,
    #                                       'nigeria_cross_imm_sus_med': 0.30000000000000004,
    #                                       'nigeria_cross_imm_sus_high': 0.7, 'nigeria_cross_imm_sev_med': 0.7,
    #                                       'nigeria_cross_imm_sev_high': 0.8500000000000001,
    #                                       'nigeria_hi5_transform_prob': 4e-10, 'nigeria_hi5_sev_fn_k': 0.08,
    #                                       'nigeria_hi5_dur_episomal_par1': 2.5, 'nigeria_hi5_dur_episomal_par2': 8.5,
    #                                       'nigeria_hi5_rel_beta': 0.7999999999999999,
    #                                       'nigeria_ohr_transform_prob': 4e-10,
    #                                       'nigeria_ohr_sev_fn_k': 0.12000000000000001,
    #                                       'nigeria_ohr_dur_episomal_par1': 2.0, 'nigeria_ohr_dur_episomal_par2': 8.5,
    #                                       'nigeria_ohr_rel_beta': 1.15, 'rwanda_beta': 0.13,
    #                                       'rwanda_cross_imm_sus_med': 0.2, 'rwanda_cross_imm_sus_high': 0.45,
    #                                       'rwanda_cross_imm_sev_med': 0.35, 'rwanda_cross_imm_sev_high': 0.8,
    #                                       'rwanda_hi5_transform_prob': 4e-10, 'rwanda_hi5_sev_fn_k': 0.19,
    #                                       'rwanda_hi5_dur_episomal_par1': 2.0, 'rwanda_hi5_dur_episomal_par2': 8.5,
    #                                       'rwanda_hi5_rel_beta': 1.1, 'rwanda_ohr_transform_prob': 4e-10,
    #                                       'rwanda_ohr_sev_fn_k': 0.05, 'rwanda_ohr_dur_episomal_par1': 2.5,
    #                                       'rwanda_ohr_dur_episomal_par2': 9.5,
    #                                       'rwanda_ohr_rel_beta': 0.8999999999999999,
    #                                       'senegal_beta': 0.15000000000000002,
    #                                       'senegal_cross_imm_sus_med': 0.35000000000000003,
    #                                       'senegal_cross_imm_sus_high': 0.45, 'senegal_cross_imm_sev_med': 0.4,
    #                                       'senegal_cross_imm_sev_high': 0.6, 'senegal_hi5_transform_prob': 6e-10,
    #                                       'senegal_hi5_sev_fn_k': 0.08, 'senegal_hi5_dur_episomal_par1': 2.0,
    #                                       'senegal_hi5_dur_episomal_par2': 9.0, 'senegal_hi5_rel_beta': 0.7,
    #                                       'senegal_ohr_transform_prob': 6e-10,
    #                                       'senegal_ohr_sev_fn_k': 0.060000000000000005,
    #                                       'senegal_ohr_dur_episomal_par1': 2.5, 'senegal_ohr_dur_episomal_par2': 10.0,
    #                                       'senegal_ohr_rel_beta': 1.15, 'sierra leone_beta': 0.26,
    #                                       'sierra leone_cross_imm_sus_med': 0.6,
    #                                       'sierra leone_cross_imm_sus_high': 0.35,
    #                                       'sierra leone_cross_imm_sev_med': 0.35,
    #                                       'sierra leone_cross_imm_sev_high': 0.9,
    #                                       'sierra leone_hi5_transform_prob': 4e-10, 'sierra leone_hi5_sev_fn_k': 0.07,
    #                                       'sierra leone_hi5_dur_episomal_par1': 2.5,
    #                                       'sierra leone_hi5_dur_episomal_par2': 7.5, 'sierra leone_hi5_rel_beta': 1.25,
    #                                       'sierra leone_ohr_transform_prob': 4e-10,
    #                                       'sierra leone_ohr_sev_fn_k': 0.060000000000000005,
    #                                       'sierra leone_ohr_dur_episomal_par1': 2.0,
    #                                       'sierra leone_ohr_dur_episomal_par2': 7.0, 'sierra leone_ohr_rel_beta': 1.05,
    #                                       'somalia_beta': 0.11, 'somalia_cross_imm_sus_med': 0.4,
    #                                       'somalia_cross_imm_sus_high': 0.3, 'somalia_cross_imm_sev_med': 0.35,
    #                                       'somalia_cross_imm_sev_high': 0.55, 'somalia_hi5_transform_prob': 2e-10,
    #                                       'somalia_hi5_sev_fn_k': 0.07, 'somalia_hi5_dur_episomal_par1': 2.5,
    #                                       'somalia_hi5_dur_episomal_par2': 8.5, 'somalia_hi5_rel_beta': 0.85,
    #                                       'somalia_ohr_transform_prob': 6e-10, 'somalia_ohr_sev_fn_k': 0.11,
    #                                       'somalia_ohr_dur_episomal_par1': 2.0, 'somalia_ohr_dur_episomal_par2': 6.5,
    #                                       'somalia_ohr_rel_beta': 1.2, 'south africa_beta': 0.15000000000000002,
    #                                       'south africa_cross_imm_sus_med': 0.5,
    #                                       'south africa_cross_imm_sus_high': 0.6000000000000001,
    #                                       'south africa_cross_imm_sev_med': 0.65,
    #                                       'south africa_cross_imm_sev_high': 0.8,
    #                                       'south africa_hi5_transform_prob': 3e-10, 'south africa_hi5_sev_fn_k': 0.07,
    #                                       'south africa_hi5_dur_episomal_par1': 2.0,
    #                                       'south africa_hi5_dur_episomal_par2': 7.5, 'south africa_hi5_rel_beta': 1.2,
    #                                       'south africa_ohr_transform_prob': 4e-10, 'south africa_ohr_sev_fn_k': 0.14,
    #                                       'south africa_ohr_dur_episomal_par1': 2.0,
    #                                       'south africa_ohr_dur_episomal_par2': 9.5, 'south africa_ohr_rel_beta': 1.1,
    #                                       'south sudan_beta': 0.17, 'south sudan_cross_imm_sus_med': 0.2,
    #                                       'south sudan_cross_imm_sus_high': 0.35, 'south sudan_cross_imm_sev_med': 0.7,
    #                                       'south sudan_cross_imm_sev_high': 0.8500000000000001,
    #                                       'south sudan_hi5_transform_prob': 2e-10, 'south sudan_hi5_sev_fn_k': 0.18,
    #                                       'south sudan_hi5_dur_episomal_par1': 2.0,
    #                                       'south sudan_hi5_dur_episomal_par2': 5.0, 'south sudan_hi5_rel_beta': 1.1,
    #                                       'south sudan_ohr_transform_prob': 2e-10, 'south sudan_ohr_sev_fn_k': 0.05,
    #                                       'south sudan_ohr_dur_episomal_par1': 2.5,
    #                                       'south sudan_ohr_dur_episomal_par2': 7.0, 'south sudan_ohr_rel_beta': 1.1,
    #                                       'tanzania_beta': 0.13, 'tanzania_cross_imm_sus_med': 0.4,
    #                                       'tanzania_cross_imm_sus_high': 0.45, 'tanzania_cross_imm_sev_med': 0.45,
    #                                       'tanzania_cross_imm_sev_high': 0.8500000000000001,
    #                                       'tanzania_hi5_transform_prob': 6e-10, 'tanzania_hi5_sev_fn_k': 0.05,
    #                                       'tanzania_hi5_dur_episomal_par1': 2.5, 'tanzania_hi5_dur_episomal_par2': 4.0,
    #                                       'tanzania_hi5_rel_beta': 0.75, 'tanzania_ohr_transform_prob': 2e-10,
    #                                       'tanzania_ohr_sev_fn_k': 0.11, 'tanzania_ohr_dur_episomal_par1': 2.5,
    #                                       'tanzania_ohr_dur_episomal_par2': 9.0, 'tanzania_ohr_rel_beta': 0.85,
    #                                       'togo_beta': 0.25, 'togo_cross_imm_sus_med': 0.30000000000000004,
    #                                       'togo_cross_imm_sus_high': 0.7, 'togo_cross_imm_sev_med': 0.45,
    #                                       'togo_cross_imm_sev_high': 0.75, 'togo_hi5_transform_prob': 4e-10,
    #                                       'togo_hi5_sev_fn_k': 0.05, 'togo_hi5_dur_episomal_par1': 3.0,
    #                                       'togo_hi5_dur_episomal_par2': 6.5, 'togo_hi5_rel_beta': 0.75,
    #                                       'togo_ohr_transform_prob': 2e-10, 'togo_ohr_sev_fn_k': 0.07,
    #                                       'togo_ohr_dur_episomal_par1': 2.0, 'togo_ohr_dur_episomal_par2': 5.0,
    #                                       'togo_ohr_rel_beta': 0.85, 'uganda_beta': 0.15000000000000002,
    #                                       'uganda_cross_imm_sus_med': 0.35000000000000003,
    #                                       'uganda_cross_imm_sus_high': 0.65, 'uganda_cross_imm_sev_med': 0.3,
    #                                       'uganda_cross_imm_sev_high': 0.5, 'uganda_hi5_transform_prob': 5e-10,
    #                                       'uganda_hi5_sev_fn_k': 0.060000000000000005,
    #                                       'uganda_hi5_dur_episomal_par1': 2.5, 'uganda_hi5_dur_episomal_par2': 6.5,
    #                                       'uganda_hi5_rel_beta': 0.8999999999999999, 'uganda_ohr_transform_prob': 3e-10,
    #                                       'uganda_ohr_sev_fn_k': 0.08, 'uganda_ohr_dur_episomal_par1': 3.0,
    #                                       'uganda_ohr_dur_episomal_par2': 10.0, 'uganda_ohr_rel_beta': 0.95,
    #                                       'zambia_beta': 0.3, 'zambia_cross_imm_sus_med': 0.45,
    #                                       'zambia_cross_imm_sus_high': 0.6000000000000001,
    #                                       'zambia_cross_imm_sev_med': 0.5, 'zambia_cross_imm_sev_high': 0.9,
    #                                       'zambia_hi5_transform_prob': 2e-10, 'zambia_hi5_sev_fn_k': 0.05,
    #                                       'zambia_hi5_dur_episomal_par1': 2.5, 'zambia_hi5_dur_episomal_par2': 8.5,
    #                                       'zambia_hi5_rel_beta': 1.2, 'zambia_ohr_transform_prob': 2e-10,
    #                                       'zambia_ohr_sev_fn_k': 0.05, 'zambia_ohr_dur_episomal_par1': 2.5,
    #                                       'zambia_ohr_dur_episomal_par2': 9.5, 'zambia_ohr_rel_beta': 0.75,
    #                                       'zimbabwe_beta': 0.25, 'zimbabwe_cross_imm_sus_med': 0.55,
    #                                       'zimbabwe_cross_imm_sus_high': 0.65, 'zimbabwe_cross_imm_sev_med': 0.5,
    #                                       'zimbabwe_cross_imm_sev_high': 0.5, 'zimbabwe_hi5_transform_prob': 2e-10,
    #                                       'zimbabwe_hi5_sev_fn_k': 0.16999999999999998,
    #                                       'zimbabwe_hi5_dur_episomal_par1': 3.0, 'zimbabwe_hi5_dur_episomal_par2': 10.0,
    #                                       'zimbabwe_hi5_rel_beta': 0.7999999999999999,
    #                                       'zimbabwe_ohr_transform_prob': 3e-10, 'zimbabwe_ohr_sev_fn_k': 0.13,
    #                                       'zimbabwe_ohr_dur_episomal_par1': 2.5, 'zimbabwe_ohr_dur_episomal_par2': 9.0,
    #                                       'zimbabwe_ohr_rel_beta': 0.7}
    #