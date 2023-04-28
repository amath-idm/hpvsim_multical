"""
This script produces figure 1 of the HPVsim calibration paper
"""

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns
import hpvsim as hpv
import hpvsim.utils as hpu

import utils as ut

#%% Run as a script
if __name__ == '__main__':

    # Settings
    n_samples = 100

    # Parameters
    mean = 54
    std = 13.5

    # Import data on cancer cases by age in Ethiopia
    cases = pd.read_csv('data/ethiopia_cancer_cases.csv')
    age_bins = cases['age'].unique()
    age_bins = np.concatenate([age_bins, [100]])
    n_cases = cases['value'].sum()
    n_bins = len(age_bins)

    # hist_data = np.zeros((n_bins-1, n_samples))
    # for n_sample in range(n_samples):
    #     synth_data = hpu.sample(dist='normal', par1=mean, par2=std, size=n_cases)
    #     hist, edges = np.histogram(synth_data, bins=age_bins)
    #     hist_data[:, n_sample] = hist
    #

    # Decompose
    debuts = [
        dict(dist='normal', par1=19.5, par2=1.5),
        dict(dist='normal', par1=19.5, par2=1.5),
    ]

    causals = [
        dict(dist='normal_pos', par1=15, par2=5),
        dict(dist='normal_pos', par1=5, par2=2)
    ]
    ttcs = [
        dict(dist='normal_pos', par1=20, par2=12),
        dict(dist='normal_pos', par1=30, par2=14)
    ]

    ut.set_font(18)
    fig, axes = pl.subplots(2, 4, figsize=(18,8))
    # for n_sample in range(n_samples):
    #     ax.plot(age_bins[:-1], hist_data[:, n_sample], c='b', alpha=0.2)

    colors = sc.gridcolors(4)
    for row in [0,1]:

        age_debut = hpu.sample(**debuts[row], size=n_cases)
        tt_causal = hpu.sample(**causals[row], size=n_cases)
        tt_cancer = hpu.sample(**ttcs[row], size=n_cases)
        age_cancer = age_debut + tt_causal + tt_cancer
        decomposed_hist, edges = np.histogram(age_cancer, bins=age_bins)

        axes[row,0].hist(age_debut, facecolor=colors[0], alpha=0.75)
        axes[row,0].set_title('Age of debut')
        axes[row,0].set_yticks([])
        axes[row,0].set_xlim([15,25])

        axes[row,1].hist(age_debut+tt_causal, facecolor=colors[1], alpha=0.75)
        axes[row,1].set_title('Age causal')
        axes[row,1].set_yticks([])
        axes[row,1].set_xlim([15,50])

        axes[row,2].hist(tt_cancer, facecolor=colors[2], alpha=0.75)
        axes[row,2].set_title('Infection duration')
        axes[row,2].set_yticks([])
        axes[row,2].set_xlim([0,75])

        axes[row,3].plot(cases['age'], cases['value'], c='r', marker='D', lw=3, label='Data')
        axes[row,3].plot(age_bins[:-1], decomposed_hist, c='k', marker='o', lw=3, label='Model 1')
        axes[row,3].set_title('Age of cancer')
        axes[row,3].set_yticks([])
        axes[row,3].legend()

    fig.tight_layout()
    pl.savefig(f"figures/nonid.png", dpi=100)


    print('Done.')
