"""
This script plots all the unconstrained calibrations together
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
def plot_single_cals(locations, n_results=20, filestem=None):
    ut.set_font(12)
    n_plots = len(locations)
    n_rows, n_cols = sc.get_rows_cols(n_plots)

    fig, axes = pl.subplots(n_rows, n_cols, figsize=(11, 10))
    axes = axes.flatten()
    resname = 'cancers'
    plot_count = 0
    date = 2020

    for pn, location in enumerate(locations):

        dflocation = location.replace(' ', '_')

        calib = sc.loadobj(f'results/unconstrained/{dflocation}_calib_{filestem}_reduced.obj')

        # Plot settings
        ax = axes[plot_count]

        # Pull out model results and data
        reslist = calib.analyzer_results
        target_data = calib.target_data[0]
        target_data = target_data[(target_data.name == resname)]

        # Make labels
        baseres = reslist[0]['cancers']
        age_labels = [str(int(baseres['bins'][i])) + '-' + str(int(baseres['bins'][i + 1])) for i in
                      range(len(baseres['bins']) - 1)]
        age_labels.append(str(int(baseres['bins'][-1])) + '+')

        # Pull out results to plot - Legacy lines commented out
        # plot_indices = calib.df.iloc[0:n_results, 0].values
        # res = [reslist[i] for i in plot_indices]
        res = reslist

        # Plot data
        x = np.arange(len(age_labels))
        ydata = np.array(target_data.value)
        ax.scatter(x, ydata, color='k', marker='s', label='Data')

        # Construct a dataframe with things in the most logical order for plotting
        bins = []
        values = []
        for run_num, run in enumerate(res):
            bins += x.tolist()
            values += list(run[resname][date])
        modeldf = pd.DataFrame({'bins': bins, 'values': values})
        sns.boxplot(ax=ax, x='bins', y='values', data=modeldf, color='b', boxprops=dict(alpha=.4))

        # Set title and labels
        title_country = location.title()
        if title_country == 'Drc':
            title_country = 'DRC'
        if title_country == 'Cote Divoire':
            title_country = "Cote d'Ivoire"
        ax.set_title(title_country)
        ax.set_ylabel('')
        ax.set_xlabel('')
        # ax.legend()
        if pn in [0, 5, 10, 15, 20, 25]:
            ax.set_ylabel('# cancers')
        if pn in [25, 26, 27, 28, 29]:
            stride = np.arange(0, len(baseres['bins']), 2)
            ax.set_xticks(x[stride], baseres['bins'].astype(int)[stride])
        else:
            ax.set_xticks(x, [])

        plot_count += 1

    fig.tight_layout()
    sc.savefig(f"figures/SMs/fig_fit_sc.png", dpi=100)


# %% Run as a script
if __name__ == '__main__':

    filestem = 'nov06'
    locations = loc.locations

    # Complete this step if you've rerun the calibrations
    # This step takes the huge calibration files and reduces them to something small enough
    # to be easily loaded and saved to the repo
    do_shrink = False
    if do_shrink:

        for pn, location in enumerate(locations):

            dflocation = location.replace(' ', '_')
            calib = sc.loadobj(f'results/unconstrained/{dflocation}_calib_{filestem}.obj')
            cal = ut.shrink_calib(calib, n_results=50)
            sc.saveobj(f'results/unconstrained/{dflocation}_calib_{filestem}_reduced.obj', cal)

    plot_single_cals(locations, n_results=50, filestem=filestem)

    print('Done.')
