from matplotlib import pyplot as plt
import json
from pandas.io.json import json_normalize
import pandas as pd
# import re

plt.rcParams['font.sans-serif'] = 'Arial'
# input correction parameters for all investigated SNe
galactic_extinction = {'SN2018hmx': {'U': 0.206, 'B': 0.173, 'V': 0.131, 'R': 0.103, 'I': 0.072,
                       'u': 0.202, 'g': 0.157, 'r': 0.109, 'i': 0.081, 'z': 0.060,
                       'J': 0.034, 'H': 0.021, 'K': 0.014, 'L': 0.007,
                       'G': 0, 'o': 0},
                       'SN1999em': {'U': 0.176, 'B': 0.147, 'V': 0.111, 'R': 0.088, 'I': 0.061,
                       'u': 0.172, 'g': 0.134, 'r': 0.093, 'i': 0.069, 'z': 0.051,
                       'J': 0.029, 'H': 0.018, 'K': 0.012, 'L': 0.006,
                       'G': 0, 'o': 0}}
distance_modulus = {'SN2018hmx': 36.06, 'SN1999em': 30.03}
z = {'SN2018hmx': 0.037, 'SN1999em': 0.0024}
discovery_date = {'SN2018hmx': '2018-10-18 00:00:00', 'SN1999em': '1999-10-29 10:33:00'}

# make sure discovery dates to datetime format
for SN in discovery_date.keys():
    discovery_date[SN] = pd.to_datetime(discovery_date[SN])

# TODO: when I combine the two analysis files to an organized project with seperate function files, make this
# TODO 'SN data dict' something that is standard and can be passed between functions and files



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

# sort the light curve data of each SN in a dict, which also has the plotting guides for each light curve
lightcurves = {'SN2018hmx': {
                    'Gaia': {'df': gaia_phot, 'marker': 'D', 'Linestyle': 'None'},
                    'ATLAS': {'df': atlas_phot, 'marker': 's', 'Linestyle': 'None'},
                    'ZTF': {'df': ztf_phot_new, 'marker': '^', 'Linestyle': 'None'},
                    'LCO': {'df': lco_phot, 'marker': 'o', 'Linestyle': 'None'}},
                'SN1999em': {
                    'Leonard': {'df': sn1999em_phot, 'marker': 'None', 'Linestyle': '--'}}}


def make_SN_dict(SN_name, lightcurves_dict, z_dict, discovery_date_dict, distance_modulus_dict,
                 galactic_extinction_dict):
    SN_dict = {'lightcurve': lightcurves_dict[SN_name],
               'name': SN_name,
               'z': z_dict[SN_name],
               'discovery_date': discovery_date_dict[SN_name],
               'distance_modulus': distance_modulus_dict[SN_name],
               'galactic_extinction': galactic_extinction_dict[SN_name]}
    return SN_dict


SN2018hmx = make_SN_dict('SN2018hmx', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)
SN1999em = make_SN_dict('SN1999em', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)

# define colormap for plotting, the colors each filter will be presented in
colormap = {'i': 'firebrick', 'r': 'tomato', 'V': 'limegreen', 'g': 'turquoise', 'B': 'blue', 'U': 'darkorchid',
            'G': 'teal', 'o': 'orange', 'R': 'tomato', 'I': 'firebrick'}

# TODO seperate to function that gets t_from_discovery and another function that normalized to rest frame by z
# TODO after that, do the t_from_discovery only for 2018hmx but rest frame normalization for both
def add_rest_frame_days_from_discovery(SN_dict):
    discovery_date = SN_dict['discovery_date']
    z = SN_dict['z']
    lightcurve = SN_dict['lightcurve']
    for source in lightcurve.keys():
        # TODO convert the timedelta to float so the plotting doesn't round it down to days
        timedelta = (lightcurve[source]['df']['datetime'] - discovery_date) / (1 + z)
        lightcurve[source]['df']['t_from_discovery'] = timedelta.dt.days
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
def lightcurve_plot(SN_dict_list, main_SN):
    fig, ax = plt.subplots(1, figsize=(9, 6))
    ax2 = ax.twinx()
    for SN in SN_dict_list:
        lightcurve = SN['lightcurve']
        for source in lightcurve.keys():
            df = lightcurve[source]['df']
            marker = lightcurve[source]['marker']
            linestyle = lightcurve[source]['Linestyle']
            for filter in df['filter'].unique():
                if SN['name'] == main_SN:
                    distance_modulus = SN['distance_modulus']
                    df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='mag', yerr='dmag', ax=ax,
                                                        marker=marker, linestyle=linestyle,
                                                        color=colormap[filter],
                                                        markeredgecolor='k', markeredgewidth=0.3)
                df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='abs_mag', ax=ax2,
                                                    marker=marker, linestyle=linestyle,
                                                    color=colormap[filter],
                                                    markeredgecolor='k', markeredgewidth =0.3,
                                                    label=source + ' (' + filter + ')',)
    ax.set_title('Light-curve over time - 2018hmx vs 1999em', fontsize=16)
    ax.set_xlabel('Time since discovery (rest-frame days)', size=16)
    ax.set_ylabel('Apparent Magnitude', size=16)
    ax.set_ylim(14, 25)
    ax2.set_ylabel('Absolute Magnitude', size=16)
    # TODO remember that the distance module difference between the y axes is hardcoded here -
    # TODO need to find way to me this automatic
    ax2.set_ylim(14 - distance_modulus, 25 - distance_modulus)
    ax2.legend(ncol=2)
    ax.invert_yaxis()
    ax.get_legend().remove()
    ax2.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    fig.savefig('light-curve over time - 2018hmx vs 1999em' + '.png')

SN2018hmx = add_rest_frame_days_from_discovery(SN2018hmx)

for SN in [SN2018hmx, SN1999em]:
    SN = remove_glactic_extinction(SN)
    SN = add_absolute_magnitude(SN)

lightcurve_plot([SN2018hmx, SN1999em], main_SN='SN2018hmx')



# show
plt.show()
