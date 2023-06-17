""" Plot age pyramids """

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns

# Imports from this repository
import locations as loc


#%% Functions
def plot_pops(locations, percentages=True):

    n_locations = len(locations)
    n_rows, n_cols = sc.get_rows_cols(n_locations)
    fig, axes = pl.subplots(n_rows, n_cols, figsize=(8, 11))
    if n_locations > 1:
        axes = axes.flatten()
    else:
        axes = axes

    m_color = '#4682b4'
    f_color = '#ee7989'
    xlabel = 'Share of population by sex' if percentages else 'Population by sex'

    for c, location in enumerate(locations):
        dflocation = location.replace(' ', '_')
        sim = sc.loadobj(f'results/sims/{dflocation}.sim')
        a = sim.get_analyzer('age_pyramid')
        pyramid = sc.odict(a.age_pyramids)[0]
        labels = list(reversed(a.age_labels))

        bins = pyramid['bins']
        ax = axes[c]

        # Prepare data
        pydf = pd.DataFrame(pyramid)
        if percentages:
            pydf['m'] = pydf['m'] / sum(pydf['m'])
            pydf['f'] = pydf['f'] / sum(pydf['f'])
        pydf['f'] = -pydf['f']  # Reverse values for females to get on same axis

        # Start making plot
        sns.barplot(x='m', y='bins', data=pydf, order=np.flip(bins), orient='h', ax=ax, color=m_color)
        sns.barplot(x='f', y='bins', data=pydf, order=np.flip(bins), orient='h', ax=ax, color=f_color)

        datadf = a.data  # [a.data.year == float(date)]
        datadf.columns = datadf.columns.str[0]
        datadf.columns = datadf.columns.str.lower()
        if percentages:
            datadf = datadf.assign(m=datadf['m'] / sum(datadf['m']), f=datadf['f'] / sum(datadf['f']))
        datadf = datadf.assign(f=-datadf['f'])
        sns.pointplot(x='m', y='a', data=datadf, order=np.flip(bins), orient='h', ax=ax, color='k', linestyles='')
        sns.pointplot(x='f', y='a', data=datadf, order=np.flip(bins), orient='h', ax=ax, color='k', linestyles='')

        ax.set_xlabel(xlabel)
        ax.set_ylabel('')
        if c in [0, 5, 10, 15, 20, 25]:
            ax.set_yticklabels(labels[1:])
        else:
            ax.set_yticklabels([])
        ax.set_xlim([-0.4, 0.4])
        xticks = ax.get_xticks()
        if percentages:
            xlabels = [f'{abs(i):.1f}' for i in xticks]
        else:
            xlabels = [f'{sc.sigfig(abs(i), sigfigs=2, SI=True)}' for i in xticks]
        if c > 24:
            ax.set_xticks(xticks, xlabels)
        else:
            ax.set_xticks(xticks, [])
        ax.set_xlabel('')
        title_country = location.title()
        if title_country == 'Drc':
            title_country = 'DRC'
        if title_country == 'Cote Divoire':
            title_country = "Cote d'Ivoire"
        ax.set_title(title_country)

    fig.tight_layout()
    sc.savefig(f'figures/0_SMs/fig_age_pyramids.png', dpi=100)


#%% Run as a script
if __name__ == '__main__':

    locations = [ll for ll in loc.locations if ll != 'congo']
    plot_pops(locations)

    print('Done.')
