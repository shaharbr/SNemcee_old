from matplotlib import pyplot as plt
import json
from pandas.io.json import json_normalize
import pandas as pd
import re

# TODO: correct for galactic extinction
plt.rcParams['font.sans-serif'] = 'Arial'

# set discovery dates and convert to datetime format
discovery_date = pd.to_datetime('17-10-18')
z = 0.037
distance_modulus = 36.06


# import LCO photometry dataset (TSV file) for 2018hmx
lco_phot = pd.read_csv('lco_photometry.txt', sep='\t', parse_dates={'datetime': ['dateobs', 'ut']})
# import Gaia photometry dataset (CSV file) for 2018hmx
gaia_phot = pd.read_csv('gaia18doi.txt', parse_dates=['#Date'], header=1)
# import ATLAS photometry dataset (TSV file) for 2018hmx
atlas_phot = pd.read_csv('ATLAS1_ACAM1.txt', sep='\t', parse_dates=['Obs-date'])
# import ZTF photometry dataset (JSON file) for 2018hmx
file = open('ZTF18accnoli_2.json')
file = json.load(file)
ztf_phot = json_normalize(file['results'])[['candidate.wall_time', 'candidate.magpsf', 'candidate.sigmapsf', 'candidate.filter']]
ztf_phot['candidate.wall_time'] = pd.to_datetime(ztf_phot['candidate.wall_time'])


# import re-analyzed ZTF photometry dataset (ASCII file) for 2018hmx
ztf_phot_new = pd.read_csv('ztf_dr1_lightcurve.txt', sep=r'\s+', header=50, usecols=['mjd|', 'mag|', 'expid|', 'catflags|'])
ztf_phot_new.drop([0, 1, 2], inplace=True)



# standardize column names and fill missing fields
gaia_phot.rename(columns={'#Date':'datetime', 'averagemag.':'mag'}, inplace=True)
atlas_phot.rename(columns={'Obs-date':'datetime', 'Mag. / Flux':'mag', 'Err':'dmag'}, inplace=True)
ztf_phot.rename(columns={'candidate.wall_time':'datetime', 'candidate.magpsf':'mag', 'candidate.sigmapsf':'dmag', 'candidate.filter':'filter'}, inplace=True)
ztf_phot_new.rename(columns={'expid|':'datetime', 'mjd|':'mag', 'mag|':'dmag', 'catflags|':'filter'}, inplace=True)

ztf_phot_new['datetime'], ztf_phot_new['mag'], ztf_phot_new['dmag'], ztf_phot_new['filter'] = \
    ztf_phot_new['datetime'].astype('float64'), ztf_phot_new['mag'].astype('float64'),\
    ztf_phot_new['dmag'].astype('float64'), ztf_phot_new['filter'].astype('str'),
ztf_phot_new['datetime'] = pd.to_datetime(ztf_phot_new['datetime'], unit='D', origin='julian')




# fill missing columns
gaia_phot['dmag'] = 0
# sort filter names
gaia_phot['filter'] = 'G'
atlas_phot['filter'] = 'o'
ztf_phot_new['filter'].replace({'zg': 'g', 'zr':'r'}, inplace=True)
lco_phot['filter'].replace({'ip': 'i', 'gp':'g', 'rp':'r'}, inplace=True)



# plotting guides
colormap = {'i': 'firebrick', 'r': 'tomato', 'V': 'limegreen', 'g': 'turquoise', 'B': 'blue', 'U': 'darkorchid', 'G': 'teal', 'o': 'orange'}
namemap = {'Gaia': {'marker': 'D', 'df': gaia_phot}, 'ATLAS': {'marker': 's', 'df': atlas_phot}, 'ZTF': {'marker': '^', 'df': ztf_phot_new}, 'LCO': {'marker': 'o', 'df': lco_phot}}

# define function for plotting light curve

def add_rest_frame_days_from_discovery(SN_dict, discovery_date, z):
    sources = SN_dict.keys()
    for source in sources:
        timedelta = (SN_dict[source]['df']['datetime'] - discovery_date) / (1 + z)
        SN_dict[source]['df']['t_from_discovery'] = timedelta.dt.days
    return SN_dict

namemap = add_rest_frame_days_from_discovery(namemap, discovery_date, z)

def add_absolute_magnitude(SN_dict, distance_modulus):
    sources = SN_dict.keys()
    for source in sources:
        SN_dict[source]['df']['abs_mag'] = SN_dict[source]['df']['mag'] - distance_modulus
    return SN_dict

namemap = add_absolute_magnitude(namemap, distance_modulus)

def lightcurve_plot(SN_dict):
    fig, ax = plt.subplots(1, figsize=(9, 6))
    ax2 = ax.twinx()
    for name in namemap:
        df = SN_dict[name]['df']
        for filter in df['filter'].unique():
            df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='mag', yerr='dmag', marker=namemap[name]['marker'],
                                                ax=ax, linestyle='None', label=name + ' (' + filter + ')'
                                                , color=colormap[filter],  markeredgecolor='k', markeredgewidth=0.3)
            df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='abs_mag', marker=namemap[name]['marker'],
                                                ax=ax2, linestyle='None', color=colormap[filter],
                                                markeredgecolor='k', markeredgewidth =0.3)
    ax.set_title('SN2018hmx light-curve over time - overlay of all sources', fontsize=16)
    ax.set_xlabel('Time since discovery (rest-frame days)', size=16)
    ax.set_ylabel('Apparent Magnitude', size=16)
    ax.set_ylim(15, 23)
    ax2.set_ylabel('Absolute Magnitude', size=16)
    ax2.set_ylim(15 - distance_modulus, 23 - distance_modulus)
    ax.legend()
    ax.invert_yaxis()
    ax2.get_legend().remove()
    ax2.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    fig.savefig(re.sub(' | - ', '_', 'SN2018hmx light-curve over time - overlay of all sources' + '.png'))


lightcurve_plot(namemap)


    # plot light curves for each source separately (in single_figs) and in an overlay figure (overlay_fig)


# single_figs = dict(zip(namemap.keys(), [0,0,0,0,0]))
# single_ax = dict(zip(namemap.keys(), [0,0,0,0,0]))
# for name in namemap:
#     lightcurve_plot(name, overlay_ax)
    # single_figs[name], single_ax[name] = plt.subplots(1, figsize=(12, 8))
    # invert magnitude y axes for display
    # plt.gca().invert_yaxis()
    # lightcurve_plot(name, single_ax[name])
    # single_ax[name].set_title('SN2018hmx light-curve over time - '+name+' only', fontsize=16)
    # single_ax[name].set_xlabel('Date', size=16)
    # single_ax[name].set_ylabel('Magnitude', size=16)
    # single_figs[name].savefig(re.sub(' |-', '_', 'SN2018hmx light-curve over time - '+name+' only'+'.png'))


# show
plt.show()
