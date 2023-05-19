"""
This script produces figure 2 of the HPVsim calibration paper
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
import settings as set
import utils as ut


#%% Plotting functions
def plot_fig3(locations):

    ut.set_font(12)
    n_plots = len(locations)
    n_rows, n_cols = sc.get_rows_cols(n_plots)

    fig, axes = pl.subplots(n_rows, n_cols, figsize=(11,10))
    axes = axes.flatten()
    resname = 'hpv_prevalence_by_age'
    plot_count = 0
    date = 2020

    for pn, location in enumerate(locations):

        # Plot settings
        ax = axes[plot_count]
        dflocation = location.replace(' ', '_')
        res = sc.loadobj(f'results/4_msims/{dflocation}.obj')

        # Plot
        x = np.arange(len(res[resname][:,-1])) # Placeholder
        ax.plot(x, res[resname][:,-1], color='k')
        ax.fill_between(x, res.hpv_prevalence_by_age.low[:,-1], res.hpv_prevalence_by_age.high[:,-1], color='b', alpha=0.2)
        ax.set_title(f'{location.capitalize()}')
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks(x, [])
        plot_count += 1

    fig.tight_layout()
    pl.savefig(f"figures/fig3.png", dpi=100)



#%% Run as a script
if __name__ == '__main__':

    locations = ['tanzania', 'angola'] #set.locations
    plot_fig3(locations)

    print('Done.')
