"""
Explore immunocompromise vs HIV prevalence.

Requires xlrd
"""

import numpy as np
import pandas as pd
import sciris as sc
import pylab as pl

T = sc.timer()
sc.options(dpi=150)

# Settings
save  = True
animate = True
start = 1990 # Earliest year data available
end   = 2020
years = [str(y) for y in sc.inclusiverange(start, end, dtype=int)]
exclude = [  # Optionally exclude outliers
    'South Africa',
    'Zimbabwe',
    'Nigeria',
    'Madagascar',
]
use_exclude = True

# Paths
relsev_fn = sc.thispath().parent / 'results' / 'rel_sevs.obj'
hiv_fn = sc.thispath().parent / 'data' / 'API_SH.DYN.AIDS.ZS_DS2_en_excel_v2_5994969.xls' # From https://data.worldbank.org/indicator/SH.DYN.AIDS.ZS

# Handle relsev
relsev = sc.load(relsev_fn)
g = relsev.groupby('Country')
gm = g.mean()

# Combine into dataframe
rs = sc.dataframe()
rs['country'] = gm.index.values
rs['relsev'] = gm.values[:,0]

# Handle HIV
hiv = sc.dataframe.read_excel(hiv_fn, header=3)
hiv = hiv.rename(columns={'Country Name':'country', 'Country Code':'code'})
keep = ['country', 'code'] + years
hiv = hiv[keep]
if use_exclude:
    hiv = hiv[~hiv.country.isin(exclude)]

# Merge
df = pd.merge(rs, hiv, on='country')

# Plot
best = -np.inf
bestyear = np.nan
fig, (ax0, ax1) = pl.subplots(2,1, figsize=(8,8))
for year in years:
    
    # Handle top plot
    pl.sca(ax0)
    pl.cla()
    x = df[year]
    inds = sc.getvalidinds(x)
    x = x[inds]
    y = df.relsev[inds]
    m = df.country[inds]
    out = sc.linregress(x, y, full=True)
    est = out.m*x + out.b
    resids = abs(y - est)
    c = sc.vectocolor(resids, cmap='turbo')
    r2 = out.corr**2
    xmin = 0
    xmax = df.max()[3:].max()
    ymin = df.relsev.min()
    ymax = df.relsev.max()
    xlim = sc.cat(xmin, xmax)
    ylim = sc.cat(ymin, ymax)
    if r2 > best:
        best = r2
        bestyear = year
    pl.plot(xlim, out.m*xlim+out.b, 'k')
    pl.scatter(x, y, c=c, s=100, alpha=0.5)
    for ix, iy, im, ic in zip(x, y, m, c):
        if not np.isnan(ix):
            pl.text(ix, iy, im, fontdict=dict(color=ic))
    pl.xlabel('HIV prevalence')
    pl.ylabel('Relative severity')
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.title(f'{year}, R²={r2:0.2f}, best={bestyear} (R²={best:0.2f})')
    if animate:
        pl.pause(0.2) # Animate
    
    # Handle bottom plot
    pl.sca(ax1)
    pl.scatter(int(year), r2, c='k')
    pl.xlabel('Year of HIV prevalence data')
    pl.ylabel('R²')
    pl.xlim([start-1, end+1])
    
    sc.figlayout()
    

if save:
    sc.savefig('relsev-hiv-scatter.png')
T.toc()