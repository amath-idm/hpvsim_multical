"""
Plot implied natural history.
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


# %% Functions

def plot_nh(sim=None):
    # Make sims
    genotypes = ['hpv16', 'hpv18']  # , 'hi5', 'ohr']
    dur_episomal = sc.autolist()
    cancer_fns = sc.autolist()
    dur_precins = sc.autolist()
    for gi, genotype in enumerate(genotypes):
        dur_precins += sim['genotype_pars'][genotype]['dur_precin']
        cancer_fns += sim['genotype_pars'][genotype]['cancer_fn']

    ####################
    # Make figure, set fonts and colors
    ####################
    ut.set_font(size=16)
    colors = sc.gridcolors(len(genotypes))
    fig, axes = pl.subplots(1, 2, figsize=(11, 5))
    axes = axes.flatten()

    ####################
    # Make plots
    ####################
    dt = 0.25
    max_x = 30
    x = np.arange(dt, max_x + dt, dt)
    annual_x = np.arange(1, 11, 1)
    width = 0.4  # the width of the bars
    multiplier = 0

    # Panel A: clearance rates
    for gi, genotype in enumerate(genotypes):
        offset = width * multiplier
        sigma, scale = ut.lognorm_params(dur_precins[gi]['par1'], dur_precins[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        axes[0].bar(annual_x + offset - width / 2, rv.cdf(annual_x), color=colors[gi], lw=2, label=genotype[-2:],
                    width=width)
        multiplier += 1
    axes[0].set_title("Clearance probabilities")
    axes[0].set_xticks(annual_x)
    axes[0].set_ylabel("Probability")
    axes[0].set_xlabel("Years since infection")
    axes[0].set_ylim([0,1])
    axes[0].legend()

    # Panel B: transform prob
    res = sim.analyzers[1].results
    years = sim.analyzers[1].durations

    df = pd.DataFrame()
    persisted = res["total"] - res["cleared"] - res["dead_other"]
    df["years"] = years
    df["prob_persisted"] = (res["persisted"]+res["progressed"]) / persisted
    df["prob_cancer"] = (res["cancer"]+res["dead_cancer"]) / persisted
    prob_cancer = (res["cancer"]+res["dead_cancer"]) / persisted

    def moving_average(a, n=3):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    # axes[1].plot(years[1:-1], moving_average(prob_cancer), lw=2, color='k')
    # axes[1].set_title("Share of all infections\n progressed to cancer")
    # axes[1].set_ylabel("Probability")
    # axes[1].set_ylim([0,1])

    ind = 31
    bottom = np.zeros(len(df["years"][:ind]))
    layers = ["prob_persisted", "prob_cancer"]
    labels = ["Persistant infection", "Cancer"]
    for ln, layer in enumerate(layers):
        axes[1].fill_between(
            df["years"][:ind],
            bottom,
            bottom + df[layer][:ind],
            label=labels[ln],
        )
        bottom += df[layer][:ind]
    axes[1].legend(loc="lower left")
    axes[1].set_title('Cancer probabilities')
    axes[1].set_xlabel("Years since infection")

    # for gi, genotype in enumerate(genotypes):
    #     cancer = hppar.compute_severity(x, pars=cancer_fns[gi])
    #     axes[1].plot(x, cancer, color=colors[gi], lw=2, label=genotype.upper())
    # axes[1].set_title("Probability of cancer\n within X years")
    # axes[1].set_ylabel("Probability")
    # axes[1].set_xlabel("Years")
    # axes[1].legend()

    # # Panel C: total dwelltime
    # dd = pd.DataFrame()
    # dw = sc.autolist()
    # gen = sc.autolist()
    # for gi, genotype in enumerate(genotypes):
    #     a = sim.get_analyzer('dwelltime_by_genotype')
    #     dw += a.dwelltime['total'][gi]
    #     gen += [genotype.upper()] * len(a.dwelltime['total'][gi])
    # dd['genotype'] = gen
    # dd['dwelltime'] = dw
    # sns.violinplot(data=dd, x="genotype", y="dwelltime", ax=axes[2], palette=colors)
    # axes[2].set_xlabel('')
    # axes[2].set_ylabel('')
    # axes[2].set_ylabel("Years")
    # axes[2].set_title('Total dwelltime\n from infection to cancer')
    #
    fig.tight_layout()

    pl.savefig(f"figures/fig3.png", dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    location = 'angola' #.replace(' ', '_')
    make_sim = False
    if make_sim:
        sim = rs.run_sim(location, ressubfolder='unconstrained', calib_par_stem=f'_pars_nov06_un',
                         analyzers=[ut.dwelltime_by_genotype(), ut.outcomes_by_year(start_year=2000)], age_pyr=True, verbose=0.1, do_save=True)
    else:
        sim = sc.loadobj(f"results/{location.replace(' ', '_')}.sim")

    plot_nh(sim)

    print('Done.')
