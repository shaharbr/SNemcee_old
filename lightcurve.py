from matplotlib import pyplot as plt
import json
from pandas.io.json import json_normalize
import pandas as pd
import re

plt.rcParams['font.sans-serif'] = 'Arial'
# TODO find real galactic extinction for 1999em (these are just made up)
galactic_extinction = {'SN2018hmx': {'U': 0.206, 'B': 0.173, 'V': 0.131, 'R': 0.103, 'I': 0.072,
                       'u': 0.202, 'g': 0.157, 'r': 0.109, 'i': 0.081, 'z': 0.060,
                       'J': 0.034, 'H': 0.021, 'K': 0.014, 'L': 0.007,
                       'G': 0, 'o': 0},
                       'SN1999em': {'U': -6.206, 'B': -6.173, 'V': -6.131, 'R': -6.103, 'I': -6.072,
                       'u': -6.202, 'g': -6.157, 'r': -6.109, 'i': -6.081, 'z': -6.060,
                       'J': -6.034, 'H': -6.021, 'K': -6.014, 'L': -6.007,
                       'G': 0, 'o': 0}}

# initiate SN dictionaries
# TODO: when I combine the two analysis files to an organized project with seperate function files, make this
# TODO 'SN data dict' something that is standard and can be passed between functions and files
# remember that some of the SN1999em data were copied from SN2018hm as placeholders, need to find its real data
SN2018hmx = {'lightcurve': {}, 'name': '2018hmx', 'z': 0.037, 'discovery_date': '2018-10-18 00:00:00',
             'distance_modulus': 36.06, 'galactic_extinction': galactic_extinction['SN2018hmx']}
SN1999em= {'lightcurve': {}, 'name': '1999em', 'z': 0, 'discovery_date': '1999-10-29 10:33:00',
           'distance_modulus': 36.06, 'galactic_extinction': galactic_extinction['SN1999em']}

# convert discovery date to datetime format
for SN in [SN2018hmx, SN1999em]:
    SN['discovery_date'] = pd.to_datetime(SN['discovery_date'])



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
# import photometry dataset (TSV file) for 1999em
sn1999em_phot = pd.read_csv('sn1999em_UBVRI_leonard02.txt', sep='\s+', header=None,
                       names=['t_from_discovery', 'U', 'B', 'V', 'R', 'I'], parse_dates=['t_from_discovery'])
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

# restructure sn1999em data as the others - all mag in one column, and another column which indicates which filter its from
sn1999em_phot = sn1999em_phot.melt(id_vars=['t_from_discovery'], var_name='filter', value_name='mag')
sn1999em_phot['t_from_discovery'] = sn1999em_phot['t_from_discovery'].astype('float64')
sn1999em_phot = sn1999em_phot.loc[sn1999em_phot['t_from_discovery'] < 300]
sn1999em_phot = sn1999em_phot.loc[sn1999em_phot['mag'] < 90]
sn1999em_phot['dmag'] = 0

# fill missing columns
gaia_phot['dmag'] = 0
# sort filter names
gaia_phot['filter'] = 'G'
atlas_phot['filter'] = 'o'
ztf_phot_new['filter'].replace({'zg': 'g', 'zr':'r'}, inplace=True)
lco_phot['filter'].replace({'ip': 'i', 'gp':'g', 'rp':'r'}, inplace=True)

# integrate the dataframes in the SN discts and add plotting guidelines
SN2018hmx['lightcurve'] = {'Gaia': {'marker': 'D', 'Linestyle': 'None', 'df': gaia_phot}, 'ATLAS': {'marker': 's', 'Linestyle': 'None', 'df': atlas_phot}, 'ZTF': {'marker': '^', 'Linestyle': 'None', 'df': ztf_phot_new}, 'LCO': {'marker': 'o', 'Linestyle': 'None', 'df': lco_phot}}
SN1999em['lightcurve'] = {'Leonard': {'marker': 'None', 'Linestyle': '--', 'df': sn1999em_phot}}

colormap = {'i': 'firebrick', 'r': 'tomato', 'V': 'limegreen', 'g': 'turquoise', 'B': 'blue', 'U': 'darkorchid',
            'G': 'teal', 'o': 'orange', 'R': 'tomato', 'I': 'firebrick'}


def add_rest_frame_days_from_discovery(SN_dict):
    discovery_date = SN_dict['discovery_date']
    z = SN_dict['z']
    lightcurve = SN_dict['lightcurve']
    for source in lightcurve.keys():
        # TODO convert the timedelta to float so the plotting doesn't round it down to days
        timedelta = (lightcurve[source]['df']['datetime'] - discovery_date) / (1 + z)
        lightcurve[source]['df']['t_from_discovery'] = timedelta.dt.days
        print(timedelta)
    SN_dict['lightcurve'] = lightcurve
    return SN_dict


def remove_glactic_extinction(SN_dict):
    lightcurve = SN_dict['lightcurve']

    galactic_extinction_values = SN_dict['galactic_extinction']
    for source in lightcurve.keys():
        df = lightcurve[source]['df']
        for filter in df['filter'].unique():
            df.loc[df['filter'] == filter, 'mag'] = df.loc[df['filter'] == filter, 'mag'] - galactic_extinction_values[filter]
    SN_dict['lightcurve'] = lightcurve
    return SN_dict


def add_absolute_magnitude(SN_dict):
    lightcurve = SN_dict['lightcurve']
    distance_modulus = SN_dict['distance_modulus']
    for source in lightcurve.keys():
        lightcurve[source]['df']['abs_mag'] = lightcurve[source]['df']['mag'] - distance_modulus
    SN_dict['lightcurve'] = lightcurve
    return SN_dict

# colormap
def lightcurve_plot(SN_dict_list):
    fig, ax = plt.subplots(1, figsize=(9, 6))
    ax2 = ax.twinx()
    for SN in SN_dict_list:
        lightcurve = SN['lightcurve']
        for source in lightcurve.keys():
            df = lightcurve[source]['df']
            marker = lightcurve[source]['marker']
            linestyle = lightcurve[source]['Linestyle']
            for filter in df['filter'].unique():
                df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='mag', yerr='dmag',
                                                    marker=marker,
                                                    ax=ax, linestyle=linestyle, label=source + ' (' + filter + ')',
                                                    color=colormap[filter],  markeredgecolor='k', markeredgewidth=0.3)
                df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='abs_mag',
                                                    marker=marker,
                                                    ax=ax2, linestyle=linestyle, color=colormap[filter],
                                                    markeredgecolor='k', markeredgewidth =0.3)
    ax.set_title('light-curve over time - overlay of all sources', fontsize=16)
    ax.set_xlabel('Time since discovery (rest-frame days)', size=16)
    ax.set_ylabel('Apparent Magnitude', size=16)
    ax.set_ylim(15, 25)
    ax2.set_ylabel('Absolute Magnitude', size=16)
    # TODO remember that the distance module difference between the y axes is hardcoded here -
    # TODO need to find way to me this automatic
    ax2.set_ylim(15 - 36.06, 25 - 36.06)
    ax.legend(ncol=2, loc='lower right')
    ax.invert_yaxis()
    ax2.get_legend().remove()
    ax2.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    # fig.savefig(re.sub(' | - ', '_', 'SN2018hmx light-curve over time - overlay of all sources' + '.png'))

SN2018hmx = add_rest_frame_days_from_discovery(SN2018hmx)

for SN in [SN2018hmx, SN1999em]:
    SN = remove_glactic_extinction(SN)
    SN = add_absolute_magnitude(SN)

lightcurve_plot([SN2018hmx, SN1999em])



# show
plt.show()
