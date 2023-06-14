''' Plot posterior distributions '''

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns

# Imports from this repository
import locations as set
import utils as ut

#%% Functions
def make_posterior_df(locations=None, n_results=50):
    dfs = sc.autolist()
    for location in locations:
        dflocation = location.replace(' ', '_')
        calib = sc.loadobj(f'results/{dflocation}_calib_apr28.obj')
        df = sc.dcp(calib.df[:n_results])
        df['location'] = location
        dfs += df
    alldf = pd.concat(dfs)
    sc.saveobj(f'results/calib_dfs_sc.obj', alldf)
    return alldf


def plot_calib_pars(locations=None):

    to_plot = [
        'hpv16_transform_prob',
        'hpv18_transform_prob',
        'hrhpv_transform_prob',
        'hrhpv_sev_fn_k',
    ]
    labels = ['Transform prob']*3 + ['k']
    alldf = sc.loadobj(f'results/calib_dfs_sc.obj')
    plotdf = alldf.loc[(alldf["location"].isin(locations))]

    # Create figure
    ut.set_font(16)
    n_rows, n_cols = sc.get_rows_cols(len(to_plot))
    fig, axes = pl.subplots(n_rows, n_cols, figsize=(8, 8))
    axes = axes.flatten()

    for pn,par in enumerate(to_plot):
        ax = axes[pn]
        legendon = True if pn==0 else False

        ax = sns.kdeplot(data=plotdf, x=par, hue="location", ax=ax, legend=False)
        # if pn==0:
        #     ax.legend(frameon=False).set_title(None)
        ax.set_title(par[:5].upper(), y=1.05)
        ax.set_xlabel(labels[pn])
        ax.set_ylabel("")

    fig.tight_layout()
    pl.savefig(f'figures/0_parplot.png', dpi=100)

    return

def plot_scatter(locations=None):

    to_plot = [
        'hpv16_transform_prob',
        'hpv18_transform_prob',
    ]
    alldf = sc.loadobj(f'results/calib_dfs_sc.obj')
    plotdf = alldf.loc[(alldf["location"].isin(locations))]

    # Create figure
    ut.set_font(16)
    fig, ax = pl.subplots(1, 1, figsize=(8, 8))
    sns.scatterplot(data=plotdf, x="hpv16_transform_prob", y="hpv18_transform_prob", hue="location", ax=ax)
    fig.tight_layout()
    pl.savefig(f'figures/0_scatter.png', dpi=100)


def plot_stripplot(locations=None):

    to_plot = [
        'hpv16_transform_prob',
        'hpv18_transform_prob',
    ]
    alldf = sc.loadobj(f'results/calib_dfs_sc.obj')
    plotdf = alldf.loc[(alldf["location"].isin(locations))]
    pdf = pd.melt(plotdf, id_vars=['index', 'mismatch', 'location'],
                  value_vars=['hpv16_transform_prob', 'hpv18_transform_prob'],
                  var_name='par', value_name='value')

    # Create figure
    ut.set_font(16)

    sns.catplot(data=pdf, x="location", y="value", row="par", color=sc.gridcolors(10), jitter=False, sharey=False, height=4, aspect=3,)
    pl.xticks(rotation=90)
    pl.savefig(f'figures/0_stripplot.png', dpi=100)


#%% Run as a script
if __name__ == '__main__':

    locations = set.locations
    # alldf = make_posterior_df(locations, n_results=10)
    # plot_calib_pars(locations)
    # plot_scatter(locations)
    plot_stripplot(locations)

    print('Done.')
