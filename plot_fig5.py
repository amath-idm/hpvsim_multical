"""
Compare natural histories from unconstrained calibrations
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


#%% Plotting function


def plot_calib_comp(locations, calib_pars=None):

    # Make sims
    genotypes = ['hpv16', 'hpv18', 'hi5', 'ohr']
    # sims = sc.objdict()
    dur_episomal = sc.autolist()
    transform_probs = sc.autolist()
    sev_fns = sc.autolist()
    dur_precins = sc.autolist()
    sims = sc.loadobj('results/comp_sims_may08.obj')
    for li, location in enumerate(locations):
        sim = sims[location]
        # pars = sc.dcp(calib_pars[li])
        # sim = hpv.Sim(pars, location=location, genotypes=genotypes)
        # sim.initialize()
        dur_precins += sim['genotype_pars']['hpv16']['dur_precin']
        dur_episomal += sim['genotype_pars']['hpv16']['dur_episomal']
        transform_probs += sim['genotype_pars']['hpv16']['transform_prob']
        sev_fns += sim['genotype_pars']['hpv16']['sev_fn']
        # sims[location] = sim

    ####################
    # Make figure, set fonts and colors
    ####################
    ut.set_font(size=16)
    colors = sc.gridcolors(len(locations))
    fig, axes = pl.subplots(1, 3, figsize=(11, 5))
    axes = axes.flatten()

    ####################
    # Make plots
    ####################
    dt = 0.25
    max_x = 30
    x = np.arange(dt, max_x+dt, dt)
    n_samples = 10
    annual_x = np.arange(1, 10, 1)
    width = 0.35  # the width of the bars
    multiplier = 0

    # Panel A: clearance rates
    for li, location in enumerate(locations):
        offset = width * multiplier
        sigma, scale = ut.lognorm_params(dur_episomal[li]['par1'], dur_episomal[li]['par2'])
        rv = lognorm(sigma, 0, scale)
        axes[0].bar(annual_x+offset-width/2, rv.cdf(annual_x), color=colors[li], lw=2, label=location.capitalize(), width=width)
        multiplier += 1
    axes[0].set_title("Proportion of clearance\n within X years")
    axes[0].set_ylabel("Probability")
    axes[0].set_xlabel("Years")
    axes[0].set_xticks(annual_x)

    # axes[0].legend()

    # Panel B: transform prob
    for li, location in enumerate(locations):
        cum_dysp = np.zeros_like(x)
        # dysp = np.zeros_like(x)
        # dysp_start_ind = sc.findnearest(x, dur_precins[li]['par1'])
        # if dysp_start_ind > 0:
        #     dysp[dysp_start_ind:] = hppar.compute_severity(x[:-dysp_start_ind], pars=sev_fns[li])
        # else:
        #     dysp[:] = hppar.compute_severity(x[:], pars=sev_fns[li])
        # axes[3].plot(x, dysp, color=colors[li], lw=2, label=location.capitalize())
        # axes[3].set_title('Dysplasia')

        dysp_int = hppar.compute_severity_integral(x, pars=sev_fns[li])
        dysp_start_ind = sc.findnearest(x, dur_precins[li]['par1'])
        if dysp_start_ind > 0:
            cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
        else:
            cum_dysp[:] = dysp_int[:]
        tp_array = hpu.transform_prob(transform_probs[li], cum_dysp)
        axes[1].plot(x, tp_array, color=colors[li], lw=2, label=location.capitalize())
    axes[1].set_title("Probability of cancer\n within X years")
    axes[1].set_ylabel("Probability")
    axes[1].set_xlabel("Years")
    axes[1].legend()

    # Panel C: total dwelltime
    dd = pd.DataFrame()
    dw = sc.autolist()
    loc = sc.autolist()
    for li, location in enumerate(locations):
        a = sims[location].get_analyzer('dwelltime_by_genotype')
        dw += a.dwelltime['total'][0]
        loc += [location.capitalize()] * len(a.dwelltime['total'][0])
    dd['location'] = loc
    dd['dwelltime'] = dw
    sns.violinplot(data=dd, x="location", y="dwelltime", ax=axes[2], cut=0, palette=colors)
    axes[2].set_xlabel('')
    axes[2].set_ylabel('')
    axes[2].set_ylabel("Years")
    axes[2].set_title('Total dwelltime\n from infection to cancer')

    fig.tight_layout()

    pl.savefig(f"figures/fig4.png", dpi=100)

    return calib_pars, sims
 

#%% Run as a script
if __name__ == '__main__':

    locations = ['tanzania', 'uganda']
    calib_pars_0 = sc.loadobj(f'results/3_sc/{locations[0]}_pars_may08_sc.obj')
    calib_pars_1 = sc.loadobj(f'results/3_sc/{locations[1]}_pars_may08_sc.obj')
    calib_pars = [calib_pars_0, calib_pars_1]
    calib_pars, sim = plot_calib_comp(locations=locations, calib_pars=calib_pars)

    print('Done.')
