from matplotlib import pyplot as plt
import json
from pandas.io.json import json_normalize
import pandas as pd
import re

print('hello')

# import LCO dataset (TSV file)
lco_phot = pd.read_csv('lco_photometry.txt', sep='\t', parse_dates={'datetime': ['dateobs', 'ut']})
# import Gaia dataset (CSV file)
gaia_phot = pd.read_csv('gaia18doi.txt', parse_dates=['#Date'], header=1)
# import ATLAS dataset (TSV file)
atlas_phot = pd.read_csv('ATLAS1_ACAM1.txt', sep='\t', parse_dates=['Obs-date'])
# import ZTF dataset (JSON file)
file = open('ZTF18accnoli_2.json')
file = json.load(file)
ztf_phot = json_normalize(file['results'])[['candidate.wall_time', 'candidate.magpsf', 'candidate.sigmapsf', 'candidate.filter']]
ztf_phot['candidate.wall_time'] = pd.to_datetime(ztf_phot['candidate.wall_time'])

# standardize column names and fill missing fields
gaia_phot.rename(columns={'#Date':'datetime', 'averagemag.':'mag'}, inplace=True)
atlas_phot.rename(columns={'Obs-date':'datetime', 'Mag. / Flux':'mag', 'Err':'dmag'}, inplace=True)
ztf_phot.rename(columns={'candidate.wall_time':'datetime', 'candidate.magpsf':'mag', 'candidate.sigmapsf':'dmag', 'candidate.filter':'filter'}, inplace=True)

# fill missing columns
gaia_phot['dmag'] = 0
gaia_phot['filter'] = 'G'
atlas_phot['filter'] = 'orange'

# plotting guides
colormap = {'ip': 'darkred', 'rp': 'orange', 'V': 'green', 'gp': 'cyan', 'B': 'blue', 'U': 'purple', 'r': 'red', 'g': 'lime', 'G': 'limegreen', 'orange': 'sandybrown'}
namemap = {'LCO': {'marker': '.', 'df': lco_phot}, 'Gaia': {'marker': '*', 'df': gaia_phot}, 'ATLAS': {'marker': 's', 'df': atlas_phot}, 'ZTF': {'marker': '^', 'df': ztf_phot}}

# define function for plotting light curve
def lightcurve_plot(name, ax):
    df = namemap[name]['df']
    for filter in df['filter'].unique():
        df.loc[df['filter'] == filter].plot(x='datetime', y='mag', yerr='dmag', marker=namemap[name]['marker'], ax=ax,
                                            linestyle='None', label=name + '_' + filter, color=colormap[filter], markersize=10)
    ax.set_title('SN2018hmx light-curve over time - overlay of all sources', fontsize=16)
    ax.set_xlabel('Date', size=16)
    ax.set_ylabel('Magnitude', size=16)

    ax.legend()


# plot light curves for each source separately (in single_figs) and in an overlay figure (overlay_fig)
overlay_fig, overlay_ax = plt.subplots(1, figsize=(12, 8))
plt.gca().invert_yaxis()

single_figs = dict(zip(namemap.keys(), [0,0,0,0]))
single_ax = dict(zip(namemap.keys(), [0,0,0,0]))
for name in namemap:
    lightcurve_plot(name, overlay_ax)
    single_figs[name], single_ax[name] = plt.subplots(1, figsize=(12, 8))
    # invert magnitude y axes for display
    plt.gca().invert_yaxis()
    lightcurve_plot(name, single_ax[name])
    single_ax[name].set_title('SN2018hmx light-curve over time - '+name+' only', fontsize=16)
    single_ax[name].set_xlabel('Date', size=16)
    single_ax[name].set_ylabel('Magnitude', size=16)
    single_figs[name].savefig(re.sub(' |-', '_', 'SN2018hmx light-curve over time - '+name+' only'+'.png'))

overlay_fig.savefig(re.sub(' | - ', '_','SN2018hmx light-curve over time - overlay of all sources'+'.png'))

# show
plt.show()
