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



#%% Plotting function

def lognorm_params(par1, par2):
    """
    Given the mean and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    mean = np.log(par1 ** 2 / np.sqrt(par2 ** 2 + par1 ** 2))  # Computes the mean of the underlying normal distribution
    sigma = np.sqrt(np.log(par2 ** 2 / par1 ** 2 + 1))  # Computes sigma for the underlying normal distribution

    scale = np.exp(mean)
    shape = sigma
    return shape, scale


def plot_sevs(calib_pars=None, genotype='hpv16'):

    cp = calib_pars['genotype_pars'][genotype]
    rel_sevs = [0.5, 0.75, 1, 1.25, 1.5]

    ####################
    # Make figure, set fonts and colors
    ####################
    ut.set_font(size=25)
    rscolors = sc.vectocolor(len(rel_sevs))
    maincolor = sc.gridcolors(1)[0]
    fig, axes = pl.subplots(2, 2, figsize=(16, 16))

    dt = 0.25
    max_x = 30
    thisx = np.arange(dt, max_x+dt, dt)

    # Panel A
    sigma, scale = lognorm_params(cp['dur_episomal']['par1'], cp['dur_episomal']['par2'])
    rv = lognorm(sigma, 0, scale)
    axes[0,0].plot(thisx, rv.pdf(thisx), color=maincolor, lw=2, label=genotype)

    # Panel C: dysplasia
    for rs,rel_sev in enumerate(rel_sevs):

        dysp = np.zeros_like(thisx)
        dysp_start_ind = sc.findnearest(thisx, cp['dur_precin']['par1'])
        if dysp_start_ind > 0:
            dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], rel_sev=rel_sev, pars=cp['sev_fn'])
        else:
            dysp[:] = hppar.compute_severity(thisx[:], rel_sev=rel_sev, pars=cp['sev_fn'])
        axes[0,1].plot(thisx, dysp, color=rscolors[rs], lw=1, label=f'{rel_sev=}')

        # Panel B: cumulative dysplasia
        cum_dysp = np.zeros_like(thisx)
        dysp_int = hppar.compute_severity_integral(thisx, rel_sev=rel_sev, pars=cp['sev_fn'])
        if dysp_start_ind > 0:
            cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
        else:
            cum_dysp[:] = dysp_int[:]
        axes[1,0].plot(thisx, dysp, color=rscolors[rs], lw=1)

        # Panel D: transform probability
        tp_array = hpu.transform_prob(cp['transform_prob'], cum_dysp)
        axes[1,1].plot(thisx, tp_array, color=rscolors[rs], lw=2, label=genotype)


    axes[0,0].set_ylabel("")
    axes[0,0].grid()
    axes[0,0].set_xlabel("Duration of infection (years)")
    axes[0,0].set_ylabel("Density")

    axes[0,1].set_ylabel("Severity of infection")
    axes[0,1].set_xlabel("Duration of infection (years)")
    axes[0,1].set_ylim([0,1])
    axes[0,1].grid()
    axes[0,1].legend()

    axes[1,0].set_ylabel("Cumulative severity of infection")
    axes[1,0].set_xlabel("Duration of infection (years)")
    axes[1,0].set_ylim([0,1])
    axes[1,0].grid()


    axes[1,1].grid()
    axes[1,1].set_ylabel("Probability of transformation")
    axes[1,1].set_xlabel("Duration of infection (years)")

    fig.tight_layout()
    pl.savefig(f"{ut.figfolder}/rel_sevs.png", dpi=100)

    return


#%% Run as a script
if __name__ == '__main__':

    location='india'
    filename = rs.mc_filename
    calib_pars = sc.loadobj(f'results/{location}_{filename}_pars.obj')
    plot_sevs(calib_pars=calib_pars, genotype='hpv16')

    print('Done.')
