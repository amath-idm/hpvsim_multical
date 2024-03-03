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
import locations as loc
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
            loclabel = location.title()
            if location == 'drc': loclabel='DRC'
            if location == 'cote divoire': loclabel = "Cote d'Ivoire"
            loc += [loclabel]*n_results
        df['Country'] = loc
        df['Mean immunocompromise level'] = rel_sevs
        sc.saveobj('results/rel_sevs.obj', df)
    else:
        df = sc.loadobj('results/rel_sevs.obj')

    hiv_prev = pd.read_csv('data/hiv_prev.csv')
    hiv_prev.set_index('country', inplace=True)
    df.set_index("Country", inplace=True)
    df = pd.merge(df, hiv_prev, left_index=True, right_index=True)
    df["country"] = df.index

    grp_order = df.groupby('country')['Mean immunocompromise level'].agg('mean').sort_values().index
    sns.boxplot(data=df, x="country", y="Mean immunocompromise level", ax=ax, order=grp_order)
    ax.yaxis.grid(True)
    pl.xticks(rotation=90)
    pl.legend([], [], frameon=False)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_ylim([0.4, 1.6])
    ax.annotate('More immunocompromised', xy=(5, 1.1), xytext=(5, 1.5), horizontalalignment="center", arrowprops=dict(arrowstyle='<-',lw=1))
    fig.tight_layout()
    sc.savefig(f"figures/fig4.png", dpi=100)

    # Correlation plots
    fig, ax = pl.subplots(1, 1, figsize=(8, 6))
    ut.set_font(24)

    # dfrs = df.groupby(['country']).mean()
    dfrs = df.groupby('country')['Mean immunocompromise level'].agg('mean')
    df2 = pd.merge(dfrs, hiv_prev, left_index=True, right_index=True)

    sns.regplot(data=df2, x="hiv_prev_scaled", y="Mean immunocompromise level", ax=ax)
#    ax.annotate('More immunocompromised', xy=(17, 1.5), xytext=(17, 1.75) ,horizontalalignment="center", arrowprops=dict(arrowstyle='<-',lw=1))
    ax.set_ylabel('Mean immunocompromise level')
    ax.set_xlabel('HIV prevalence among women 15-49')

    fig.tight_layout()
    sc.savefig(f"figures/SMs/fig_hiv_scatter.png", dpi=100)

    # ax = axes[1]
    # sns.regplot(data=df2, x="le", y="Mean immunocompromise level", ax=axes[1])
    # ax.set_ylabel('')
    # ax.set_xlabel('Life expectancy')


    return df2


#%% Run as a script
if __name__ == '__main__':

    locations = loc.locations
    df2 = plot_rel_sevs(locations, make=False, filestem='_calib_nov06_iv')

    print('Done.')
