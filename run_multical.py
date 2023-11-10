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
import pandas as pd
import utils as ut

# Imports from this repository
import run_sim as rs
import calibration as cal
import locations as loc

# Comment out to not run
to_run = [
    'run_calibration',
    # 'plot_calibration',
]

debug = False  # Smaller runs
do_save = True


# Run settings for calibration (dependent on debug)
n_trials    = [120, 1][debug]  # How many trials to run for calibration
n_workers   = [40, 1][debug]    # How many cores to use
storage     = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations

########################################################################
# Run calibration
########################################################################
def make_unique_priors(locations=None):
    ''' Make priors for the parameters that vary across settings '''

    unique_pars = dict()
    for location in locations:
        unique_pars[location] = dict(
            calib_pars = dict(
                beta=[0.2, 0.1, 0.34, 0.02],
                m_cross_layer=[0.3, 0.1, 0.7, 0.05],
                m_partners=dict(
                    c=dict(par1=[0.2, 0.1, 0.6, 0.02])
                ),
                f_cross_layer=[0.1, 0.05, 0.5, 0.05],
                f_partners=dict(
                    c=dict(par1=[0.2, 0.1, 0.6, 0.02])
                ),
            ),
            genotype_pars = dict(
                hi5=dict(
                    cin_fn=dict(k=[.2, .15, .4, 0.01]),
                ),
                ohr=dict(
                    cin_fn=dict(k=[.2, .15, .4, 0.01]),
                ),
            )
        )

    return unique_pars

def make_posterior_df(locations=None, n_results=50):
    dfs = sc.autolist()
    for location in locations:
        dflocation = location.replace(' ', '_')
        calib = sc.loadobj(f'results/unconstrained/{dflocation}_calib_nov06.obj')
        df = sc.dcp(calib.df[:n_results])
        df['location'] = location
        dfs += df
    alldf = pd.concat(dfs)
    sc.saveobj(f'results/calib_dfs_sc.obj', alldf)
    return alldf


def run_calib(locations=None, sc_pars=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):

    # Define shared calibration parameters - same values used across sims
    common_pars = dict(
        genotype_pars=dict(
            hpv16=dict(
                cin_fn=dict(k=[.35, .25, .45, 0.01]),
            ),
            hpv18=dict(
                cin_fn=dict(k=[.35, .25, .45, 0.01]),
            )
        )
    )

    if sc_pars is None:
        unique_pars = make_unique_priors(locations)
    else:
        unique_pars = None

    sims = []
    for location in locations:
        if sc_pars is not None:
            calib_pars = sc_pars[location]
        else:
            calib_pars = None
        sim = rs.make_sim(location, calib_pars=calib_pars)
        sim.label = location
        sims.append(sim)

    calib = cal.MultiCal(
        sims,
        common_pars=common_pars,
        unique_pars=unique_pars,
        name=f'multical1106',
        datafiles=ut.make_datafiles(locations),
        load_if_exists=True,
        db_name='multical1106.db',
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

    calib = sc.load(f'results/constrained/multical{filestem}.obj')

    if save_pars:
        sims = []
        for location in locations:
            pars_file = f'results/constrained/{location}_multical{filestem}_pars.obj'
            calib_pars = calib.trial_pars_to_sim_pars(slabel=location, which_pars=which_pars)
            sc.save(pars_file, calib_pars)

    if do_plot:
        sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
        sc.options(font='Libertinus Sans')
        for location in locations:
            fig = calib.plot(slabel=location, res_to_plot=50, plot_type='sns.boxplot')
            fig.suptitle(f'Calibration results, {location.capitalize()}')
            fig.tight_layout()
            fig.savefig(f'figures/constrained/multical{filestem}_{location}.png')


    return calib, sims


#%% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    filestem = '_nov06'
    locations = ['angola', 'benin']  # loc.locations

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        sc_pars = None
        sims, calib = run_calib(locations=locations, sc_pars=sc_pars, n_trials=n_trials, n_workers=n_workers, do_save=do_save, do_plot=False, filestem=filestem)

        # Get the parameters from the single cals
        # alldf = make_posterior_df(locations, n_results=1)
        # alldf = sc.loadobj('results/calib_dfs_sc.obj')
        # sc_pars = dict()
        # for location in locations:
        #     thisdf = alldf[alldf.location == location]
        #     sc_pars[location] = dict(
        #         genotype_pars=dict(
        #             hi5=dict(
        #                 cin_fn=dict(k=thisdf.hi5_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50),
        #                 cancer_fn=dict(ld50=thisdf.hi5_cancer_fn_ld50.iloc[0], method='cin_integral', **dict(k=thisdf.hi5_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50)),
        #                 dur_cin=dict(dist='lognormal', par1=thisdf.hi5_dur_cin_par1.iloc[0], par2=thisdf.hi5_dur_cin_par2.iloc[0]),
        #                 rel_beta=thisdf.hi5_rel_beta.iloc[0]
        #             ),
        #             ohr=dict(
        #                 cin_fn=dict(k=thisdf.ohr_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50),
        #                 cancer_fn=dict(ld50=thisdf.ohr_cancer_fn_ld50.iloc[0], method='cin_integral', **dict(k=thisdf.ohr_cin_fn_k.iloc[0], form='logf2', x_infl=0, ttc=50)),
        #                 dur_cin=dict(dist='lognormal', par1=thisdf.ohr_dur_cin_par1.iloc[0], par2=thisdf.ohr_dur_cin_par2.iloc[0]),
        #                 rel_beta=thisdf.ohr_rel_beta.iloc[0]
        #             )
        #         )
        #     )
        #     if ~np.isnan(thisdf.cross_imm_sev_high.iloc[0]):
        #         sc_pars[location]['cross_imm_sev_high'] = thisdf.cross_imm_sev_high.iloc[0]
        #         sc_pars[location]['cross_imm_sev_med'] = thisdf.cross_imm_sev_med.iloc[0]
        #         sc_pars[location]['cross_imm_sus_high'] = thisdf.cross_imm_sus_high.iloc[0]
        #         sc_pars[location]['cross_imm_sus_med'] = thisdf.cross_imm_sus_med.iloc[0]


    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        calib, sims = load_calib(filestem=filestem, locations=locations, do_plot=True, save_pars=True) # lo_hiv_locations

    T.toc('Done')