"""
Read in sexual behavior data
"""

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns
from scipy.stats import norm, lognorm
import hpvsim as hpv

# Imports from this repository
import utils as ut
import run_sim as rs

def percentiles_to_pars(x1, p1, x2, p2):
    """ Find the parameters of a normal distribution where:
            P(X < p1) = x1
            P(X < p2) = x2
    """
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)

    location = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = (x2 - x1) / (p2ppf - p1ppf)
    return location, scale


def logn_percentiles_to_pars(x1, p1, x2, p2):
    """ Find the parameters of a lognormal distribution where:
            P(X < p1) = x1
            P(X < p2) = x2
    """
    x1 = np.log(x1)
    x2 = np.log(x2)
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)
    s = (x2 - x1) / (p2ppf - p1ppf)
    mean = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = np.exp(mean)
    return s, scale


def read_debut_data(save_pars=True, dist_type='lognormal'):
    '''
    Read in dataframes taken from DHS and return them in a plot-friendly format,
    optionally saving the distribution parameters
    '''

    df1 = pd.read_csv('data/afs_dist.csv')
    df2 = pd.read_csv('data/afs_median.csv')

    # Deal with median data
    df2['y'] = 50

    # Rearrange data into a plot-friendly format
    dff = {}
    rvs = {'Women':[], 'Men':[]}

    for sex in ['Women', 'Men']:

        dfw = df1[['Country', f'{sex} 15', f'{sex} 18', f'{sex} 20', f'{sex} 22', f'{sex} 25', f'{sex} never']]
        dfw = dfw.melt(id_vars='Country', value_name='Percentage', var_name='AgeStr')

        # Add values for proportion ever having sex
        countries = dfw.Country.unique()
        n_countries = len(countries)
        vals = []
        for country in countries:
            val = 100-dfw.loc[(dfw['AgeStr'] == f'{sex} never') & (dfw['Country'] == country) , 'Percentage'].iloc[0]
            vals.append(val)

        data_cat = {'Country': countries, 'AgeStr': [f'{sex} 60']*n_countries}
        data_cat["Percentage"] = vals
        df_cat = pd.DataFrame.from_dict(data_cat)
        dfw = pd.concat([dfw,df_cat])

        conditions = [
            (dfw['AgeStr'] == f"{sex} 15"),
            (dfw['AgeStr'] == f"{sex} 18"),
            (dfw['AgeStr'] == f"{sex} 20"),
            (dfw['AgeStr'] == f"{sex} 22"),
            (dfw['AgeStr'] == f"{sex} 25"),
            (dfw['AgeStr'] == f"{sex} 60"),
        ]
        values = [15, 18, 20, 22, 25, 60]
        dfw['Age'] = np.select(conditions, values)

        dff[sex] = dfw

        # Optionally save the distribution parameters
        if save_pars:
            res = dict()
            res["location"] = []
            res["par1"] = []
            res["par2"] = []
            res["dist"] = []
            for pn,country in enumerate(countries):
                dfplot = dfw.loc[(dfw["Country"] == country) & (dfw["AgeStr"] != f'{sex} never') & (dfw["AgeStr"] != f'{sex} 60')]
                x1 = 15
                p1 = dfplot.loc[dfplot["Age"] == x1, 'Percentage'].iloc[0] / 100
                x2 = 25
                p2 = dfplot.loc[dfplot["Age"] == x2, 'Percentage'].iloc[0] / 100
                res["location"].append(country)
                res["dist"].append(dist_type)

                if dist_type=='normal':
                    loc, scale = percentiles_to_pars(x1, p1, x2, p2)
                    rv = norm(loc=loc, scale=scale)
                    res["par1"].append(loc)
                    res["par2"].append(scale)
                elif dist_type=='lognormal':
                    s, scale = logn_percentiles_to_pars(x1, p1, x2, p2)
                    rv = lognorm(s=s, scale=scale)
                    res["par1"].append(rv.mean())
                    res["par2"].append(rv.std())

                rvs[sex].append(rv)

            pd.DataFrame.from_dict(res).to_csv(f'data/sb_pars_{sex.lower()}.csv')

    return countries, dff, df2, rvs


def read_marriage_data():
    '''
    Read in dataframes taken from DHS and return them in a plot-friendly format,
    optionally saving the distribution parameters
    '''

    dfraw = pd.read_csv('data/prop_married.csv')
    df = dfraw.melt(id_vars=['Country', 'Survey'], value_name='Percentage', var_name='AgeRange')
    return df


def get_sb_from_sims(which='prop_married'):
    '''
    Run sims with the sexual debut parameters inferred from DHA data, and save
    the proportion of people of each age who've ever had sex
    '''

    locations = ut.locations
    dataless_locations = ut.nosbdata_locations
    data_locations = [loc for loc in locations if loc not in dataless_locations]
    countries_to_run = locations
    sims = rs.run_sims(locations=countries_to_run, analyzers=[ut.AFS(),ut.prop_married()], debug=False)

    # Prepare to save model output
    dfs = sc.autolist()

    for country in countries_to_run:
        a = sims[country].get_analyzer(which)
        if which=='AFS':
            for cs,cohort_start in enumerate(a.cohort_starts):
                df = pd.DataFrame()
                df['age'] = a.bins
                df['cohort'] = cohort_start
                df['model_prop_f'] = a.prop_active_f[cs,:]
                df['model_prop_m'] = a.prop_active_m[cs,:]
                df['country'] = country
                dfs += df
        elif which=='prop_married':
            df = a.df
            df['country'] = country
            dfs += df

    alldf = pd.concat(dfs)
    sc.saveobj(f'results/model_sb_{which}.obj', alldf)

    return alldf


def plot_sb(make_sb_pars=True, dist_type='lognormal'):
    '''
    Create plots of sexual behavior inputs and outputs
    '''

    countries, dff, df2, rvs = read_debut_data(save_pars=make_sb_pars, dist_type=dist_type)
    alldf = sc.loadobj(f'results/model_sb_AFS.obj')
    n_countries = len(countries)
    n_rows, n_cols = sc.get_rows_cols(n_countries)

    for sk,sex in {'f':'Women', 'm':'Men'}.items():
        fig, axes = pl.subplots(n_rows, n_cols, figsize=(8,11))
        axes = axes.flatten()
        dfw = dff[sex]

        for pn,country in enumerate(countries):
            ax = axes[pn]
            dfplot = dfw.loc[(dfw["Country"]==country)&(dfw["AgeStr"]!=f'{sex} never')&(dfw["AgeStr"]!=f'{sex} 60')]
            dfmed = df2.loc[df2["Country"] == country]
            sns.scatterplot(ax=ax, data=dfplot, x="Age", y="Percentage")
            sns.scatterplot(ax=ax, data=dfmed, x=f"{sex} median", y="y")

            rv = rvs[sex][pn]
            xx = np.arange(12,30.1,0.1)
            xxx = np.arange(12,31,1)
            ax.plot(xx, rv.cdf(xx)*100, 'k--', lw=2)

            location = ut.rev_map_sb_loc(country)
            for cohort in alldf["cohort"].unique():
                modely = np.array(alldf.loc[(alldf["country"]==location)&(alldf["cohort"]==cohort)][f'model_prop_{sk}'])
                ax.plot(xxx, modely*100, 'b-', lw=1)
            title_country = country
            if title_country == 'Congo Democratic Republic':
                title_country = 'DRC'

            ax.set_title(title_country)
            ax.set_ylabel('')
            ax.set_xlabel('')

        fig.tight_layout()
        pl.savefig(f"figures/fig_sb_{sex.lower()}.png", dpi=100)

    return


def plot_prop_married():
    '''
    Create plots of sexual behavior inputs and outputs
    '''

    # Read in data and model results
    df = read_marriage_data()
    modeldf = sc.loadobj(f'results/model_sb_prop_married.obj')
    modeldf.reset_index()

    # Plot
    countries = df.Country.unique()
    n_countries = len(countries)
    n_rows, n_cols = sc.get_rows_cols(n_countries)
    colors = sc.gridcolors(1)

    fig, axes = pl.subplots(n_rows, n_cols, figsize=(8, 11))
    axes = axes.flatten()

    for pn, country in enumerate(countries):
        ax = axes[pn]

        # Plot data
        dfplot_d = df.loc[(df["Country"] == country)]
        sns.scatterplot(ax=ax, data=dfplot_d, x="AgeRange", y="Percentage")

        # Plot model
        location = ut.rev_map_sb_loc(country)
        dfplot_m = modeldf.loc[modeldf['country']==location]
        # dfplot_m['val'] = dfplot_m['val'].multiply(100)
        dfplot_m['val'] = dfplot_m['val'].apply(lambda x: x * 100)
        sns.boxplot(data=dfplot_m, x="age", y="val", color=colors[0], ax=ax)

        title_country = country
        if title_country == 'Congo Democratic Republic':
            title_country = 'DRC'

        ax.set_title(title_country)
        ax.set_ylabel('')
        ax.set_xlabel('')

    fig.tight_layout()
    pl.savefig(f"figures/fig_prop_married.png", dpi=100)

    return




#%% Run as a script
if __name__ == '__main__':


    # plot_sb(make_sb_pars=True, dist_type='lognormal')
    # alldf = prop_ever_from_sims()
    df = get_sb_from_sims('prop_married')
    # plot_prop_married()

    print('Done.')
