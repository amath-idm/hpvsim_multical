"""![](figures/fig4_calib_comp.png)
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
def plot_rel_sevs(locations, make=False, filestem=None, n_results=100):

    ut.set_font(16)
    fig, ax = pl.subplots(1, 1, figsize=(11,5))

    if make:
        df = pd.DataFrame()
        rel_sevs = sc.autolist()
        loc = sc.autolist()
        for pn, location in enumerate(locations):
            dflocation = location.replace(' ', '_')
            calib = sc.loadobj(f'results/1a_iv/{dflocation+filestem}.obj')
            rel_sevs += [calib.trial_pars_to_sim_pars(which_pars=i)['sev_dist']['par1'] for i in range(n_results)]
            loclabel = location.capitalize()
            if location=='drc':loclabel='DRC'
            if location=='cote divoire': loclabel="Cote d'Ivoire"
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
    ax.annotate('More immunocompromised', xy=(5, 1.2), xytext=(5, 1.6) ,horizontalalignment="center", arrowprops=dict(arrowstyle='<-',lw=1))
    fig.tight_layout()
    pl.savefig(f"figures/fig3.png", dpi=100)



#%% Run as a script
if __name__ == '__main__':

    locations = set.locations
    plot_rel_sevs(locations, make=False, filestem='_calib_may18')

    print('Done.')
