"""
This script produces figure 1 of the HPVsim calibration paper
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


#%% Plotting functions
def plot_fig1(locations, calib, n_results=20):

    n_plots = len(locations)
    n_rows, n_cols = sc.get_rows_cols(n_plots)

    fig, axes = pl.subplots(n_rows, n_cols, figsize=(8,11))
    axes = axes.flatten()
    resname = 'cancers'
    plot_count = 0
    date = 2020

    for pn, location in enumerate(locations):

        # Plot settings
        ax = axes[plot_count]

        # Pull out model results and data
        reslist = calib.age_results[location]
        target_data = calib.target_data[location][0]
        target_data = target_data[(target_data.name == resname)]

        # Make labels
        baseres = reslist[0]['cancers']
        age_labels = [str(int(baseres['bins'][i])) + '-' + str(int(baseres['bins'][i + 1])) for i in range(len(baseres['bins'])-1)]
        age_labels.append(str(int(baseres['bins'][-1])) + '+')

        # Pull out results to plot
        plot_indices = calib.df.iloc[0:n_results, 0].values
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

    locations = rm.locations
    calib = sc.loadobj('results/multical_apr03.obj')
    plot_fig1(locations, calib, n_results=20)

    print('Done.')
