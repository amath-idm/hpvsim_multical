'''
This script explore the relationship between rel_sev and transform_prob
'''

# Import packages
import sciris as sc
import pylab as pl
import hpvsim.parameters as hppar
import hpvsim.utils as hpu
import hpvsim as hpv
import numpy as np

# Imports from this repository
import utils as ut

#%% Functions
def explore_rel_sev(location='kenya'):

    # Group genotypes
    genotypes = ['hpv16', 'hpv18', 'hrhpv']
    sim = hpv.Sim(location=location, genotypes=genotypes)
    sim.initialize()
    ng = len(genotypes)
    genotype_pars = sim['genotype_pars']

    dur_precin = [genotype_pars[g]['dur_precin'] for g in range(ng)]
    transform_probs = [genotype_pars[g]['transform_prob'] for g in range(ng)]
    sev_fns = [genotype_pars[g]['sev_fn'] for g in range(ng)]

    dt = 0.25
    max_x = 30
    thisx = np.arange(dt, max_x+dt, dt)
    gi=0
    rel_sevs = np.arange(0.5,1.6,.1)
    transform_prob_vars = np.arange(2e-11,29e-11,2.5e-11)
    colors = sc.vectocolor(len(rel_sevs))

    # Make figure, set fonts and colors
    ut.set_font(size=12)
    fig, axes = pl.subplots(1, 3, figsize=(11, 3))

    # Make outputs
    for rs,rel_sev in enumerate(rel_sevs):

        # Plot dysplasia
        dysp = np.zeros_like(thisx)
        dysp_start_ind = sc.findnearest(thisx, dur_precin[gi]['par1'])
        if dysp_start_ind > 0:
            dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], rel_sev=rel_sev, pars=sev_fns[gi])
        else:
            dysp[:] = hppar.compute_severity(thisx[:], pars=sev_fns[gi])

        axes[0].plot(thisx, dysp, color=colors[rs], lw=2, label=f'{rel_sev=}')
        axes[0].set_title('Severity curve')

        # Plot transform prob
        cum_dysp = np.zeros_like(thisx)
        dysp_int = hppar.compute_severity_integral(thisx, rel_sev=rel_sev, pars=sev_fns[gi])
        if dysp_start_ind > 0:
            cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
        else:
            cum_dysp[:] = dysp_int[:]
        tp_array = hpu.transform_prob(transform_probs[gi], cum_dysp)

        axes[1].plot(thisx, tp_array, color=colors[rs], lw=2, label=f'{rel_sev=}')
        axes[1].set_title('Transform prob')

    # Plot different transform probs
    # Dysp with rel_sev=1
    cum_dysp = np.zeros_like(thisx)
    dysp_int = hppar.compute_severity_integral(thisx, rel_sev=1, pars=sev_fns[gi])
    if dysp_start_ind > 0:
        cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
    else:
        cum_dysp[:] = dysp_int[:]

    for t,tp in enumerate(transform_prob_vars):
        tp_array = hpu.transform_prob(tp, cum_dysp)
        axes[2].plot(thisx, tp_array, color=colors[t], lw=2, label=f'{tp=}')
        axes[2].set_title('Transform prob')

    fig.tight_layout()
    pl.savefig(f"figures/0_rel_sev.png", dpi=100)


#%% Run as a script
if __name__ == '__main__':

    explore_rel_sev()