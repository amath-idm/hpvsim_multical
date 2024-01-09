"""
This script plots the type distributions
"""

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns

# Imports from this repository
import locations as loc
import utils as ut


# %% Plotting functions
def plot_types(locations, mc_calib=None, n_results=20, filestem=None, which='unconstrained'):
    ut.set_font(12)
    n_plots = len(locations)
    glabels = ['16', '18', 'Hi5', 'OHR']
    fig, axes = sc.getrowscols(n_plots, figsize=(11, 10), make=True, remove_extra=True)

    axes = axes.flatten()
    gen_cols = sc.gridcolors(4)
    ms = 80

    plot_count = 0

    for pn, location in enumerate(locations):

        dflocation = location.replace(' ', '_')

        if mc_calib is None:
            calib = sc.loadobj(f'results/{which}/{dflocation}_calib_{filestem}.obj')
            datadf = calib.target_data[-1]
            reslist = calib.sim_results
        else:
            calib = mc_calib
            datadf = calib.target_data[location][-1]
            reslist = calib.sim_results[location]

        ydata = datadf.value.values
        x = np.arange(len(ydata))

        # Plot settings
        ax = axes[plot_count]

        # CINS and cancer by genotype
        rkey = 'cancerous_genotype_dist'

        # Pull out the analyzer and sim results
        index_to_plot = calib.df.iloc[0:n_results, 0].values
        sim_results = [reslist[i] for i in index_to_plot]

        # Extract model results
        bins = []
        values = []
        for run_num, run in enumerate(sim_results):
            bins += x.tolist()
            if sc.isnumber(run[rkey]):
                values += sc.promotetolist(run[rkey])
            else:
                values += run[rkey].tolist()
        modeldf = pd.DataFrame({'bins': bins, 'values': values})

        # Plot model
        sns.boxplot(ax=ax, x='bins', y='values', data=modeldf, palette=gen_cols, showfliers=True)
        ax.scatter(x, ydata, color='k', marker='d', s=ms)

        ax.set_ylim([0, 1])
        ax.set_xticks(np.arange(4), calib.glabels)
        ax.set_ylabel('')
        ax.set_xlabel('')

        # Set title and labels
        title_country = location.title()
        if title_country == 'Drc':
            title_country = 'DRC'
        if title_country == 'Cote Divoire':
            title_country = "Cote d'Ivoire"
        ax.set_title(title_country)

        plot_count += 1

        ax.set_ylim([0, 1])
        ax.set_xticks(np.arange(4), glabels)

    fig.tight_layout()
    pl.savefig(f"figures/SMs/fig_types_{which}.png", dpi=100)

 
# %% Run as a script
if __name__ == '__main__':

    locations = loc.cancer_type_locs

    calibration_stems = dict(
        unconstrained='nov06',
        immunovarying='nov06_iv',
    )
    for which, filestem in calibration_stems.items():
        plot_types(locations, n_results=50, filestem=filestem, which=which)

    mc_calib = sc.loadobj('results/constrained/multical_nov13.obj')
    plot_types(locations, mc_calib=mc_calib, n_results=50, which='constrained')

    print('Done.')
