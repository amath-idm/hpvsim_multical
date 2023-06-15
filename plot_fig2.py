"""
Plot immuno-varying calibrations
"""

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns


# Imports from this repository
import run_sim as rs
import run_multical as rm
import locations as set
import utils as ut


#%% Plotting functions
def plot_fig1(locations, filestem=None, n_results=20):

    ut.set_font(12)
    n_plots = len(locations)
    n_rows, n_cols = sc.get_rows_cols(n_plots)

    fig, axes = pl.subplots(n_rows, n_cols, figsize=(11,10))
    axes = axes.flatten()
    resname = 'cancers'
    plot_count = 0
    date = 2020

    for pn, location in enumerate(locations):

        # Plot settings
        ax = axes[plot_count]

        dflocation = location.replace(' ', '_')
        sccalib = sc.loadobj(f'results/1a_iv/{dflocation}_calib_{filestem}.obj')
        reslist = sccalib.analyzer_results
        target_data = sccalib.target_data[0]
        target_data = target_data[(target_data.name == resname)]

        # Make labels
        baseres = reslist[0]['cancers']
        age_labels = [str(int(baseres['bins'][i])) + '-' + str(int(baseres['bins'][i + 1])) for i in range(len(baseres['bins'])-1)]
        age_labels.append(str(int(baseres['bins'][-1])) + '+')

        # Pull out results to plot
        plot_indices = sccalib.df.iloc[0:n_results, 0].values
        res = [reslist[i] for i in plot_indices]

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
        sns.boxplot(ax=ax, x='bins', y='values', data=modeldf, color='b')

        # Set title and labels
        # ax.set_xlabel('Age group')
        ax.set_title(f'{location.capitalize()}')
        ax.set_ylabel('')
        ax.set_xlabel('')
        # ax.legend()
        ax.set_xticks(x, [])
        plot_count += 1

    fig.tight_layout()
    pl.savefig(f"figures/fig1.png", dpi=100)



#%% Run as a script
if __name__ == '__main__':

    locations = set.locations
    plot_fig1(locations, filestem='may18', n_results=50)

    print('Done.')
