"""
This script produces figure 4 of the HPVsim methods paper, showing the natural history
"""
import hpvsim as hpv
import hpvsim.utils as hpu
import hpvsim.parameters as hppar
import pylab as pl
import pandas as pd
from scipy.stats import lognorm, norm
import numpy as np
import sciris as sc
import utils as ut
import seaborn as sns

import run_sim as rs


#%% Functions

def plot_nh(sim=None):

    # Make sims
    genotypes = ['hpv16', 'hpv18'] #, 'hi5', 'ohr']
    dur_episomal = sc.autolist()
    transform_probs = sc.autolist()
    sev_fns = sc.autolist()
    dur_precins = sc.autolist()
    for gi, genotype in enumerate(genotypes):
        dur_precins += sim['genotype_pars'][genotype]['dur_precin']
        dur_episomal += sim['genotype_pars'][genotype]['dur_episomal']
        transform_probs += sim['genotype_pars'][genotype]['transform_prob']
        sev_fns += sim['genotype_pars'][genotype]['sev_fn']
        # sims[location] = sim

    ####################
    # Make figure, set fonts and colors
    ####################
    ut.set_font(size=16)
    colors = sc.gridcolors(len(genotypes))
    fig, axes = pl.subplots(1, 3, figsize=(11, 5))
    axes = axes.flatten()

    ####################
    # Make plots
    ####################
    dt = 0.25
    max_x = 30
    x = np.arange(dt, max_x+dt, dt)
    annual_x = np.arange(1, 11, 1)
    width = 0.4  # the width of the bars
    multiplier = 0

    # Panel A: clearance rates
    for gi, genotype in enumerate(genotypes):
        offset = width * multiplier
        sigma, scale = ut.lognorm_params(dur_episomal[gi]['par1'], dur_episomal[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        axes[0].bar(annual_x+offset-width/2, rv.cdf(annual_x), color=colors[gi], lw=2, label=genotype.upper(), width=width)
        multiplier += 1
    axes[0].set_title("Proportion clearing\n within X years")
    axes[0].set_xticks(annual_x)
    # axes[0].legend()

    # Panel B: transform prob
    for gi, genotype in enumerate(genotypes):
        cum_dysp = np.zeros_like(x)
        dysp_int = hppar.compute_severity_integral(x, pars=sev_fns[gi])
        dysp_start_ind = sc.findnearest(x, dur_precins[gi]['par1'])
        if dysp_start_ind > 0:
            cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
        else:
            cum_dysp[:] = dysp_int[:]
        tp_array = hpu.transform_prob(transform_probs[gi], cum_dysp)
        axes[1].plot(x, tp_array, color=colors[gi], lw=2, label=genotype.upper())
    axes[1].set_title("Probability of cancer\n within X years")
    axes[1].legend()

    # Panel C: total dwelltime
    dd = pd.DataFrame()
    dw = sc.autolist()
    gen = sc.autolist()
    for gi, genotype in enumerate(genotypes):
        a = sim.get_analyzer('dwelltime_by_genotype')
        dw += a.dwelltime['total'][gi]
        gen += [genotype.upper()] * len(a.dwelltime['total'][gi])
    dd['genotype'] = gen
    dd['dwelltime'] = dw
    sns.violinplot(data=dd, x="genotype", y="dwelltime", ax=axes[2], palette=colors)
    axes[2].set_xlabel('')
    axes[2].set_ylabel('')
    axes[2].set_title('Total dwelltime\n from infection to cancer')

    fig.tight_layout()

    pl.savefig(f"figures/fig2.png", dpi=100)

    return
 

#%% Run as a script
if __name__ == '__main__':

    sim = sc.loadobj('results/mali.sim')
    plot_nh(sim)

    print('Done.')
