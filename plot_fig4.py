"""
Plot distribution of rel_sevs
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
def plot_rel_sevs(locations, make=False, filestem=None, n_results=50):

    ut.set_font(16)
    # fig, ax = pl.subplots(1, 1, figsize=(11,5))
    fig, ax = pl.subplots(1, 1, figsize=(11,5))

    if make:
        df = pd.DataFrame()
        rel_sevs = sc.autolist()
        loc = sc.autolist()
        for pn, location in enumerate(locations):
            dflocation = location.replace(' ', '_')
            calib = sc.loadobj(f'results/immunovarying/{dflocation+filestem}.obj')
            rel_sevs += [calib.trial_pars_to_sim_pars(which_pars=i)['sev_dist']['par1'] for i in range(n_results)]
            loclabel = location.capitalize()
            if location == 'drc': loclabel='DRC'
            if location == 'cote divoire': loclabel = "Cote d'Ivoire"
            loc += [loclabel]*n_results
        df['Country'] = loc
        df['Mean immunocompromise level'] = rel_sevs
        sc.saveobj('results/rel_sevs.obj', df)
    else:
        df = sc.loadobj('results/rel_sevs.obj')
    grp_order = df.groupby('Country')['Mean immunocompromise level'].agg('mean').sort_values().index
    sns.boxplot(data=df, x="Country", y="Mean immunocompromise level", ax=ax, order=grp_order)
    ax.yaxis.grid(True)
    pl.xticks(rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_ylim([0, 2.5])
    ax.annotate('More immunocompromised', xy=(5, 1.2), xytext=(5, 1.6) ,horizontalalignment="center", arrowprops=dict(arrowstyle='<-',lw=1))
    fig.tight_layout()
    pl.savefig(f"figures/fig4.png", dpi=100)

    # Correlation plots
    fig, axes = pl.subplots(1, 2, figsize=(11,5))

    hiv_prev = pd.read_csv('data/hiv_prev.csv')
    dfrs = df.groupby(['Country']).mean()
    hiv_prev.set_index('country', inplace=True)
    df2 = pd.merge(dfrs, hiv_prev, left_index=True, right_index=True)

    ax = axes[0]
    sns.regplot(data=df2, x="hiv_prev_scaled", y="Mean immunocompromise level", ax=axes[0])
#    ax.annotate('More immunocompromised', xy=(17, 1.5), xytext=(17, 1.75) ,horizontalalignment="center", arrowprops=dict(arrowstyle='<-',lw=1))
    ax.set_ylabel('Mean immunocompromise level')
    ax.set_xlabel('HIV prevalence among women 15-49')

    ax = axes[1]
    sns.regplot(data=df2, x="le", y="Mean immunocompromise level", ax=axes[1])
    ax.set_ylabel('')
    ax.set_xlabel('Life expectancy')

    fig.tight_layout()
    pl.savefig(f"figures/scatter.png", dpi=100)


#%% Run as a script
if __name__ == '__main__':

    locations = set.locations
    plot_rel_sevs(locations, make=False, filestem='_calib_jun15')

    print('Done.')
