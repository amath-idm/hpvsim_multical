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

def transform_prob(tp, dysp):
    '''
    Returns transformation probability given dysplasia
    Using formula for half an ellipsoid:
        V = 1/2 * 4/3 * pi * a*b*c
          = 2 * a*b*c
          = 2* dysp * (dysp/2)**2, assuming that b = c = 1/2 a
          = 1/2 * dysp**3
    '''
    # return 1-np.power(1-tp, (dysp*100))
    return 1-np.power(1-tp, 0.5*((dysp*100)**3))



def plot_fig4(location, calib_pars=None, old_pars=True):

    # Group genotypes
    genotypes = ['hpv16', 'hpv18', 'hrhpv']
    if calib_pars is not None:
        pars = sc.dcp(calib_pars)
        sim = hpv.Sim(pars, location=location, genotypes=genotypes)
    else:
        sim = hpv.Sim(location=location, genotypes=genotypes)
    sim.initialize()

    # # Get parameters
    ng = sim['n_genotypes']
    genotype_map = sim['genotype_map']

    # Get parameters
    if old_pars:
        clinical_cutoffs = dict(precin=0.03, cin1=0.353, cin2=0.676, cin3=0.8)
    else:
        clinical_cutoffs = dict(cin1=0.33, cin2=0.676, cin3=0.8)

    if calib_pars is not None:
        genotype_pars = sim['genotype_pars']
        if old_pars:
            dur_precin = [
                dict(dist='uniform', par1=0, par2=0.05),
                dict(dist='uniform', par1=0, par2=0.05),
                dict(dist='uniform', par1=0, par2=0.05),
            ]
        else:
            dur_precin = [genotype_pars[genotype_map[g]]['dur_precin'] for g in range(ng)]
        dur_episomal = [genotype_pars[genotype_map[g]]['dur_episomal'] for g in range(ng)]
        transform_probs = [genotype_pars[genotype_map[g]]['transform_prob'] for g in range(ng)]
        sev_fns = [genotype_pars[genotype_map[g]]['sev_fn'] for g in range(ng)]
        sev_dist = sim['sev_dist']

    else:
        # Define parameters
        dur_precin = [
            dict(dist='normal_pos', par1=0.5, par2=0.25),
            dict(dist='normal_pos', par1=0.5, par2=0.25),
            dict(dist='normal_pos', par1=0.5, par2=0.25),
        ]

        dur_episomal = [
            dict(dist='lognormal', par1=2, par2=10), # HPV16
            dict(dist='lognormal', par1=2, par2=10), # 18
            dict(dist='lognormal', par1=2, par2=10), # OHR
        ]
        transform_probs = [
            9 / 1e11,
            8 / 1e11,
            7 / 1e11,
        ]

        sev_fns = [
            dict(form='logf2', k=0.3, x_infl=0, ttc=20),
            dict(form='logf2', k=0.2, x_infl=0, ttc=25),
            dict(form='logf2', k=0.15, x_infl=0, ttc=25),
        ]

        sev_dist = dict(dist='normal_pos', par1=1.0, par2=0.15)

    ####################
    # Make figure, set fonts and colors
    ####################
    ut.set_font(size=25)
    colors = sc.gridcolors(ng)
    cmap = pl.cm.Oranges([0.25, 0.5, 0.75, 1])
    fig, ax = pl.subplot_mosaic('AB;AB;CD;CE;FG;FG', figsize=(16, 24))

    ####################
    # Panel A and C
    ####################

    glabels = ['HPV16', 'HPV18', 'HRHPV']
    ####################
    # Make plots
    ####################

    dt = 0.25
    max_x = 30
    thisx = np.arange(dt, max_x+dt, dt)
    n_samples = 10

    # Durations and severity of dysplasia
    for gi, gtype in enumerate(genotypes):

        # Panel A: durations
        sigma, scale = lognorm_params(dur_episomal[gi]['par1'], dur_episomal[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        ax['A'].plot(thisx, rv.pdf(thisx), color=colors[gi], lw=2, label=glabels[gi])

        # Panel C: dysplasia
        dysp = np.zeros_like(thisx)
        dysp_start_ind = sc.findnearest(thisx, dur_precin[gi]['par1'])
        if dysp_start_ind > 0:
            dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], pars=sev_fns[gi])
        else:
            dysp[:] = hppar.compute_severity(thisx[:], pars=sev_fns[gi])
        ax['C'].plot(thisx, dysp, color=colors[gi], lw=2, label=gtype.upper())

        # Panel B: transform prob
        cum_dysp = np.zeros_like(thisx)
        if not sev_fns[gi].get('s') or sev_fns[gi]['s'] == 1:
            dysp_int = hppar.compute_severity_integral(thisx, pars=sev_fns[gi])
        else:
            cumdysp = sim['cumdysp'][genotype_map[gi]]
            t = np.around(thisx / dt).astype(int)  # Round
            t[t > len(cumdysp) - 1] = len(cumdysp) - 1
            dysp_int = cumdysp[t]

        if dysp_start_ind > 0:
            cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
        else:
            cum_dysp[:] = dysp_int[:]

        tp_array = transform_prob(transform_probs[gi], cum_dysp)
        ax['F'].plot(thisx, tp_array, color=colors[gi], lw=2, label=gtype.upper())

        # Add rel_sev samples to B and C plots
        for smpl in range(n_samples):
            rel_sev = hpu.sample(**sev_dist)
            dysp_start = hpu.sample(**dur_precin[gi])
            rel_dysp = np.zeros_like(thisx)
            rel_cum_dysp = np.zeros_like(thisx)
            dysp_start_ind = sc.findnearest(thisx, dysp_start)

            if not sev_fns[gi].get('s') or sev_fns[gi]['s'] == 1:
                rel_dysp_int = hppar.compute_severity_integral(thisx, rel_sev=rel_sev, pars=sev_fns[gi])
            else: # Don't know how to ger rel_sev in here
                cumdysp = sim['cumdysp'][genotype_map[gi]]
                t = np.around(thisx / dt).astype(int)  # Round
                t[t > len(cumdysp) - 1] = len(cumdysp) - 1
                rel_dysp_int = cumdysp[t]

            if dysp_start_ind>0:
                rel_dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], rel_sev=rel_sev, pars=sev_fns[gi])
                rel_cum_dysp[dysp_start_ind:] = rel_dysp_int[:-dysp_start_ind]
            elif dysp_start_ind==0:
                rel_dysp[:] = hppar.compute_severity(thisx[:], rel_sev=rel_sev, pars=sev_fns[gi])
                rel_cum_dysp[:] = rel_dysp_int[:]
            ax['C'].plot(thisx, rel_dysp, color=colors[gi], lw=1, alpha=0.5)

            rel_tp_array = transform_prob(transform_probs[gi], rel_cum_dysp)
            ax['F'].plot(thisx, rel_tp_array, color=colors[gi], lw=1, alpha=0.5)

    ax['A'].set_ylabel("")
    ax['A'].grid()
    ax['A'].set_xlabel("Duration of infection (years)")
    ax['A'].set_ylabel("Density")
    ax['A'].legend()

    ax['C'].set_ylabel("Severity of infection")
    ax['C'].set_xlabel("Duration of infection (years)")
    ax['C'].set_ylim([0,1])
    ax['C'].grid()

    ax['C'].axhline(y=clinical_cutoffs['cin1'], ls=':', c='k')
    ax['C'].axhline(y=clinical_cutoffs['cin2'], ls=':', c='k')
    ax['C'].axhspan(0, clinical_cutoffs['cin1'], color=cmap[0], alpha=.4)
    ax['C'].axhspan(clinical_cutoffs['cin1'], clinical_cutoffs['cin2'], color=cmap[1], alpha=.4)
    ax['C'].axhspan(clinical_cutoffs['cin2'], 1.0, color=cmap[2], alpha=.4)
    ax['C'].text(-0.3, 0.15, 'CIN1', rotation=90)
    ax['C'].text(-0.3, 0.48, 'CIN2', rotation=90)
    ax['C'].text(-0.3, 0.8, 'CIN3', rotation=90)

    ax['F'].grid()
    ax['F'].set_ylabel("Probability of transformation")
    ax['F'].set_xlabel("Duration of infection (years)")

    ####################
    # Panel D
    ####################

    # This section calculates the overall share of outcomes for people infected with each genotype
    precinshares, cin1shares, cin2shares, cin3shares, cancershares = [], [], [], [], [] # Initialize the share of people who get dysplasia vs cancer
    ttt = []
    ttcs = []
    n_cancers = []
    dw1s = []
    dw2s = []
    dw3s = []

    # Loop over genotypes
    for gi in range(ng):

        # First, determine the outcomes for women
        sigma, scale = lognorm_params(dur_episomal[gi]['par1'], dur_episomal[gi]['par2'])
        rv = lognorm(sigma, 0, scale)

        peak_dysp = np.zeros_like(thisx)
        cum_dysp = np.zeros_like(thisx)
        dysp_start_ind = sc.findnearest(thisx, dur_precin[gi]['par1'])
        if not sev_fns[gi].get('s') or sev_fns[gi]['s'] == 1:
            dysp_int = hppar.compute_severity_integral(thisx, pars=sev_fns[gi])
        else:
            cumdysp = sim['cumdysp'][genotype_map[gi]]
            t = np.around(thisx / dt).astype(int)  # Round
            t[t > len(cumdysp) - 1] = len(cumdysp) - 1
            dysp_int = cumdysp[t]
        if dysp_start_ind>0:
            peak_dysp[dysp_start_ind:] = hppar.compute_severity(thisx[:-dysp_start_ind], pars=sev_fns[gi])
            cum_dysp[dysp_start_ind:] = dysp_int[:-dysp_start_ind]
        else:
            peak_dysp[:] = hppar.compute_severity(thisx[:], pars=sev_fns[gi])
            cum_dysp[:] = dysp_int[:]

        tp = transform_prob(transform_probs[gi], cum_dysp)

        # To start find women who advance to cancer
        cancer_inds = hpu.true(hpu.n_binomial(tp, len(thisx)))

        # For each genotype, make a histogram of the distribution of times to transformation
        ttt.append(thisx[cancer_inds])
        n_cancers.append(len(cancer_inds))

        # Also calculate dwelltimes
        rel_sev_samples = hpu.sample(**sev_dist, size=len(cancer_inds))
        dur_to_cin3 = hppar.compute_inv_severity(clinical_cutoffs['cin3'], rel_sev=rel_sev_samples, pars=sev_fns[gi])
        ttc = np.maximum(thisx[cancer_inds], dur_to_cin3)
        ttcs.append(ttc)
        # dur_transformed = ttc-thisx[cancer_inds]
        dw1s.append(hppar.compute_inv_severity(clinical_cutoffs['cin1'], rel_sev=rel_sev_samples, pars=sev_fns[gi]))
        dw2s.append(hppar.compute_inv_severity(clinical_cutoffs['cin2'], rel_sev=rel_sev_samples, pars=sev_fns[gi]))
        dw3s.append(ttc-hppar.compute_inv_severity(clinical_cutoffs['cin3'], rel_sev=rel_sev_samples, pars=sev_fns[gi]))

        # Find women who only advance to PRECIN
        if clinical_cutoffs.get('precin'):
            indprecin = sc.findinds((peak_dysp < clinical_cutoffs['precin']))[-1]
            n_precin = len(sc.findinds(peak_dysp < clinical_cutoffs['precin']))
            indcin1 = sc.findinds((peak_dysp > clinical_cutoffs['precin']) & (peak_dysp < clinical_cutoffs['cin1']))[-1]
            n_cin1 = len(sc.findinds((peak_dysp > clinical_cutoffs['precin']) & (peak_dysp < clinical_cutoffs['cin1'])))
        else:
            indprecin = sc.findinds(peak_dysp == 0)[-1]
            n_precin = len(sc.findinds(peak_dysp ==0))
            indcin1 = sc.findinds((peak_dysp < clinical_cutoffs['cin1']))[-1]
            n_cin1 = len(sc.findinds((peak_dysp < clinical_cutoffs['cin1'])))

        precin_share = rv.cdf(thisx[indprecin])
        cin1_share = rv.cdf(thisx[indcin1]) - rv.cdf(thisx[indprecin])

        # See if there are women who advance to CIN2 and get their indices if so
        if (peak_dysp > clinical_cutoffs['cin1']).any():
            n_cin2 = len(sc.findinds((peak_dysp > clinical_cutoffs['cin1']) & (peak_dysp < clinical_cutoffs['cin2'])))
            indcin2 = sc.findinds((peak_dysp > clinical_cutoffs['cin1']) & (peak_dysp < clinical_cutoffs['cin2']))[-1]
        else:
            n_cin2 = 0
            indcin2 = indcin1
        cin2_share = rv.cdf(thisx[indcin2]) - rv.cdf(thisx[indcin1])

        if (peak_dysp > clinical_cutoffs['cin2']).any():
            n_cin3 = len(sc.findinds(peak_dysp > clinical_cutoffs['cin2']))
            indcin3 = sc.findinds((peak_dysp > clinical_cutoffs['cin2']))[-1]  # Index after which people develop CIN3 (plus possibly cancer)
        else:
            n_cin3 = 0
            indcin3 = indcin2
        cin3_share = rv.cdf(thisx[indcin3]) - rv.cdf(thisx[indcin2])

        if clinical_cutoffs.get('precin'):
            n_cancer_precin = len(np.intersect1d(cancer_inds, sc.findinds((peak_dysp < clinical_cutoffs['precin']))))
            n_cancer_cin1 = len(np.intersect1d(cancer_inds, sc.findinds((peak_dysp > clinical_cutoffs['precin']) & (peak_dysp < clinical_cutoffs['cin1']))))
        else:
            n_cancer_precin= len(np.intersect1d(cancer_inds, sc.findinds(peak_dysp == 0)))
            n_cancer_cin1 = len(np.intersect1d(cancer_inds, sc.findinds((peak_dysp < clinical_cutoffs['cin1']))))
        n_cancer_cin2 = len(np.intersect1d(cancer_inds, sc.findinds((peak_dysp > clinical_cutoffs['cin1']) & (peak_dysp < clinical_cutoffs['cin2']))))
        n_cancer_cin3 = len(np.intersect1d(cancer_inds, sc.findinds((peak_dysp > clinical_cutoffs['cin2']))))

        cancer_share_of_precins = n_cancer_precin/n_precin
        cancer_share_of_cin1s = n_cancer_cin1 / n_cin1  # Share of CIN1 women who get cancer
        cancer_share_of_cin2s = n_cancer_cin2 / n_cin2  # Share of CIN2 women who get cancer
        cancer_share_of_cin3s = n_cancer_cin3 / n_cin3  # Share of CIN3 women who get cancer

        precin_share *= 1 - cancer_share_of_precins
        cin1_share *= 1 - cancer_share_of_cin1s
        cin2_share *= 1 - cancer_share_of_cin2s
        cin3_share *= 1 - cancer_share_of_cin3s
        cancer_share = 1 - (precin_share + cin1_share + cin2_share + cin3_share)

        precinshares.append(precin_share)
        cin1shares.append(cin1_share)
        cin2shares.append(cin2_share)
        cin3shares.append(cin3_share)
        cancershares.append(cancer_share)

    # Final plot
    bottom = np.zeros(ng)
    all_shares = [precinshares,
                  cin1shares,
                  cin2shares,
                  cin3shares,
                  cancershares
                  ]

    for gn, grade in enumerate(['Pre-CIN', 'CIN1', 'CIN2', 'CIN3', 'Cancer']):
        ydata = np.array(all_shares[gn])
        color = cmap[gn - 1, :] if gn > 0 else 'gray'
        ax['B'].bar(np.arange(1, ng + 1), ydata, color=color, bottom=bottom, label=grade)
        bottom = bottom + ydata

    ax['B'].set_xticks(np.arange(1,ng + 1))
    ax['B'].set_xticklabels(glabels)
    ax['B'].set_ylabel("")
    ax['B'].set_ylabel("Distribution of outcomes")
    # ax['E'].legend(bbox_to_anchor=(1.1, 1))
    handles, labels = ax['D'].get_legend_handles_labels()
    ax['B'].legend(handles, labels, frameon=True, loc='lower right')

    # Distribution of times to transformation
    bins = np.linspace(0, max_x, 11)
    ax['D'].hist(ttt, label=glabels, bins=bins, color=colors)
    ax['D'].set_xlabel("Time from infection to transformation")
    ax['D'].set_ylabel("Density")

    ax['E'].hist(ttcs, label=glabels, bins=bins, color=colors)
    ax['E'].set_xlabel("Time from infection to cancer")
    ax['E'].set_ylabel("Density")

    dw1 = np.concatenate(dw1s) #[dw1s[0]] * len(dw3s[0]) + [dw1s[1]] * len(dw3s[1]) + [dw1s[2]] * len(dw3s[2])
    dw2 = np.concatenate(dw2s) #[dw2s[0]] * len(dw3s[0]) + [dw2s[1]] * len(dw3s[1]) + [dw2s[2]] * len(dw3s[2])
    dw3 = np.concatenate(dw3s)
    gnames = ['hpv16']*len(dw3s[0]) + ['hpv18']*len(dw3s[1]) + ['hrhpv']*len(dw3s[2])
    df = pd.DataFrame(dict(dw1=dw1,dw2=dw2,dw3=dw3,genotype=gnames))
    df['id'] = df.index
    df = pd.wide_to_long(df, ["dw"], i="id", j="CIN")
    df = df.reset_index()
    sns.violinplot(data=df, x="CIN", y="dw", hue="genotype", ax=ax['G'], cut=0)
    ax['G'].set_xlabel("CIN")
    ax['G'].set_ylabel("Dwelltime")

    # fs=40
    # pl.figtext(0.02, 0.955, 'A', fontweight='bold', fontsize=fs)
    # pl.figtext(0.51, 0.955, 'C', fontweight='bold', fontsize=fs)
    # pl.figtext(0.02, 0.47, 'B', fontweight='bold', fontsize=fs)
    # pl.figtext(0.51, 0.47, 'D', fontweight='bold', fontsize=fs)
    fig.tight_layout()

    pl.savefig(f"figures/calib_nathx_{location}.png", dpi=100)

    return calib_pars, sim


#%% Run as a script
if __name__ == '__main__':

    # location='nigeria_new'
    location='ethiopia'
    # filename = rs.mc_filename
    # calib_pars = sc.loadobj(f'results/{location}_{filename}_pars.obj')

    calib_pars = sc.loadobj(f'results/ethiopia_pars_apr26_sc.obj')
    # calib_pars['genotype_pars']['hpv16']['transform_prob'] = 5e-10
    # calib_pars['genotype_pars']['hpv18']['transform_prob'] = 2e-10
    # calib_pars['genotype_pars']['hrhpv']['transform_prob'] = 3e-5
    # calib_pars['genotype_pars']['hrhpv']['sev_fn']['k'] = 0.2
    # calib_pars['genotype_pars']['hpv16']['sev_fn']['k'] = 0.35
    # calib_pars['genotype_pars']['hpv16']['dur_episomal'] = {'dist': 'lognormal', 'par1': 0.88, 'par2': 0.88}
    calib_pars, sim = plot_fig4(location=location, calib_pars=calib_pars, old_pars=False)

    print('Done.')
