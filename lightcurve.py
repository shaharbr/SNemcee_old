from matplotlib import pyplot as plt
import json
from pandas.io.json import json_normalize
import pandas as pd
# import re
import numpy as np
import data_import
import lightcurve_param_fit

# TODO split to modules:
# 1) data imports
# 2) functions that organize the data in data df and SN_dicts. maybe also include in this the functions
#  that normalize and calculate t from discovery
# 3) plotting
# 4) light curve parameter calculations


# TODO:
# standardize the SN_dict object so it is the same for both the photometry and spectroscopy analyses

plt.rcParams['font.sans-serif'] = 'Arial'
# TODO ask - what should I do regarding the G and o filters which dont have galactic extinction values in the NES?

# TODO turn the data input dicts to stardardized files that can be imported with the data (csv?)
# input correction parameters for all investigated SNe
distance_modulus = {'SN2018hmx': 36.06, 'SN1999em': 30.03, 'SN2004ej': 32.3, 'SN2012A': 29.73 , 'ASASSN14kg': 34.04}
z = {'SN2018hmx': 0.037, 'SN1999em': 0.0024, 'SN2004ej': 0.0091, 'SN2012A': 0.002291429, 'ASASSN14kg': 0.014477}
# note: SN2012A z calculated as average of multiple z values onon the sn2 catalog
discovery_date = {'SN2018hmx': '2018-10-17 15:30:14', 'SN1999em': '1999-10-29 10:33:00',
                  'SN2004ej': '2004-09-10', 'SN2012A': '2012-01-07', 'ASASSN14kg': '2014-11-17'}
galactic_extinction = {'SN2018hmx': {'U': 0.206, 'B': 0.173, 'V': 0.131, 'R': 0.103, 'I': 0.072,
                       'u': 0.202, 'g': 0.157, 'r': 0.109, 'i': 0.081, 'z': 0.060,
                       'J': 0.034, 'H': 0.021, 'K': 0.014, 'L': 0.007,
                       'G': 0, 'o': 0},
                       'SN1999em': {'U': 0.176, 'B': 0.147, 'V': 0.111, 'R': 0.088, 'I': 0.061,
                       'u': 0.172, 'g': 0.134, 'r': 0.093, 'i': 0.069, 'z': 0.051,
                       'J': 0.029, 'H': 0.018, 'K': 0.012, 'L': 0.006,
                       'G': 0, 'o': 0},
                       'SN2004ej': {'V': 0},
                       'SN2012A': {'B': 0, 'V': 0, 'U': 0, 'UVW1': 0, 'UVW2': 0, 'UVM2': 0},
                       'ASASSN14kg': {'B': 0, 'V': 0, 'U': 0, 'i': 0, 'r': 0, 'g': 0}}
# define colormap for plotting, the colors each filter will be presented in
colormap = {'i': 'firebrick', 'r': 'tomato', 'g': 'turquoise',
            'V': 'limegreen', 'B': 'blue', 'U': 'darkorchid', 'G': 'teal', 'R': 'tomato', 'I': 'firebrick',
            'o': 'orange',
            'UVW1': 'darkorchid', 'UVW2': 'darkorchid', 'UVM2': 'darkorchid'}

# convert discovery dates to MJD
for SN in discovery_date.keys():
    discovery_date[SN] = data_import.convert_to_mjd(discovery_date[SN], from_datetime=True)


# import photometry data files
lco_phot = data_import.lco_phot('lco_photometry.txt')
gaia_phot = data_import.gaia_phot('gaia18doi.txt')
atlas_phot = data_import.atlas_phot('ATLAS1_ACAM1.txt')
ztf_phot_new = data_import.ztf_phot_new('ztf_dr1_lightcurve.txt')
sn1999em_leonard_phot = data_import.leonard_phot('sn1999em_UBVRI_leonard02.txt')
sn2004ej_phot = data_import.sne_catalog_phot('SN2004ej.csv')
sn2012A_phot = data_import.sne_catalog_phot('SN2012A.csv')
ASASSN14kg_phot = data_import.sne_catalog_phot('ASASSN-14kg.csv')




# sort the light curve data of each SN in a dict, which also has the plotting guides for each light curve
lightcurves = {'SN2018hmx': {
                    'Gaia': {'df': gaia_phot, 'marker': 'D', 'Linestyle': 'None'},
                    'ATLAS': {'df': atlas_phot, 'marker': 's', 'Linestyle': 'None'},
                    'ZTF': {'df': ztf_phot_new, 'marker': '^', 'Linestyle': 'None'},
                    'LCO': {'df': lco_phot, 'marker': 'o', 'Linestyle': 'None'}},
                'SN1999em': {
                    'Leonard': {'df': sn1999em_leonard_phot, 'marker': 'None', 'Linestyle': '--'}},
                'SN2004ej': {
                    'Planck': {'df': sn2004ej_phot, 'marker': 'None', 'Linestyle': '--'}},
                'SN2012A': {
                    'Swift': {'df': sn2012A_phot, 'marker': 'None', 'Linestyle': '--'}},
                'ASASSN14kg': {
                    'ASASSN': {'df': ASASSN14kg_phot, 'marker': 'None', 'Linestyle': '--'}},
                }


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
SN2004ej = make_SN_dict('SN2004ej', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)
SN2012A = make_SN_dict('SN2012A', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)
ASASSN14kg = make_SN_dict('ASASSN14kg', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)


# TODO:
# seperate to function that gets t_from_discovery and another function that normalized to rest frame by z
# after that, do the t_from_discovery only for 2018hmx but rest frame normalization for both

def add_rest_frame_days_from_discovery(SN_dict):
    discovery_date = SN_dict['discovery_date']
    z = SN_dict['z']
    lightcurve = SN_dict['lightcurve']
    for source in lightcurve.keys():
        timedelta = (lightcurve[source]['df']['mjd'] - discovery_date) / (1 + z)
        lightcurve[source]['df']['t_from_discovery'] = timedelta
    SN_dict['lightcurve'] = lightcurve
    return SN_dict

# remove before-discovery data points
def remove_data_before_discovery(SN_dict):
    lightcurve = SN_dict['lightcurve']
    for source in lightcurve.keys():
        lightcurve[source]['df'] = lightcurve[source]['df'].loc[lightcurve[source]['df']['t_from_discovery'] >= 0]
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
def lightcurve_plot(SN_dict_list, main_SN, V50_line=False):
    fig, ax = plt.subplots(1, figsize=(9, 6))
    ax2 = ax.twinx()
    for SN in SN_dict_list:
        lightcurve = SN['lightcurve']
        for source in lightcurve.keys():
            df = lightcurve[source]['df']
            marker = lightcurve[source]['marker']
            linestyle = lightcurve[source]['Linestyle']
            for filter in df['filter'].unique():
                df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='abs_mag', ax=ax2,
                                                    marker=marker, linestyle=linestyle,
                                                    color=colormap[filter],
                                                    markeredgecolor='k', markeredgewidth =0.3,
                                                    label=source + ' (' + filter + ')',)
                if SN['name'] == main_SN:
                    distance_modulus = SN['distance_modulus']
                    df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='mag', yerr='dmag', ax=ax,
                                                        marker=marker, linestyle=linestyle,
                                                        color=colormap[filter],
                                                        markeredgecolor='k', markeredgewidth=0.3)
                    if isinstance(V50_line, np.ndarray) and source == 'LCO':
                        V50_poly1d = np.poly1d(V50_line)
                        # TODO make drawing borders same 50 days as ones used for extracting fit
                        x = list(df.loc[(df['filter'] == 'V') &
                                        (df['t_from_discovery'] <= 80), 't_from_discovery'])
                        y = V50_poly1d(x)
                        ax2.plot(x, y, color='k')

    ax.set_title('Light-curve over time - '+str([SN['name'] for SN in SN_dict_list]), fontsize=16)
    ax.set_xlabel('Time since discovery (rest-frame days)', size=16)
    ax.set_ylabel('Apparent Magnitude', size=16)
    ax.set_ylim(14, 27)
    ax2.set_ylabel('Absolute Magnitude', size=16)
    # TODO remember that the distance module difference between the y axes is hardcoded here -
    # TODO need to find way to me this automatic
    ax2.set_ylim(14 - distance_modulus, 27 - distance_modulus)
    ax2.legend(ncol=2)
    ax.invert_yaxis()
    ax.get_legend().remove()
    ax2.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)

    # fig.savefig('light-curve over time - 2018hmx vs 1999em' + '.png')

SN2018hmx = add_rest_frame_days_from_discovery(SN2018hmx)
SN2004ej = add_rest_frame_days_from_discovery(SN2004ej)
SN2012A = add_rest_frame_days_from_discovery(SN2012A)
ASASSN14kg = add_rest_frame_days_from_discovery(ASASSN14kg)


for SN in [SN2018hmx, SN1999em, SN2012A, SN2004ej, ASASSN14kg]:
    SN = remove_glactic_extinction(SN)
    SN = add_absolute_magnitude(SN)

SN2018hmx = remove_data_before_discovery(SN2018hmx)
SN1999em = remove_data_before_discovery(SN1999em)


def remove_LCO_outlier(SN_dict):
    df = SN_dict['lightcurve']['LCO']['df']
    df = df.loc[df['t_from_discovery'] < 190]
    return df

SN2018hmx['lightcurve']['LCO']['df'] = remove_LCO_outlier(SN2018hmx)


lightcurve_plot([SN2018hmx, SN1999em], main_SN='SN2018hmx')
lightcurve_plot([SN2018hmx, SN2004ej], main_SN='SN2018hmx')
lightcurve_plot([SN2018hmx, SN2012A], main_SN='SN2018hmx')
lightcurve_plot([SN2018hmx, ASASSN14kg], main_SN='SN2018hmx')

# light curve parameters:

s50V_2018hmx = lightcurve_param_fit.calc_s50V(SN2018hmx)
p0_2018hmx = lightcurve_param_fit.calc_p0(SN2018hmx, time_range=[140, 190])

# lightcurve_plot([SN2018hmx], main_SN='SN2018hmx', V50_line=v50_regression_params)
print('s50V:', s50V_2018hmx)
print('p0:', p0_2018hmx)

# lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 'p0')
# lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 's50V')

# sampler = lightcurve_param_fit.SN_lightcurve_params(SN2018hmx)
# lightcurve_param_fit.chain_plots(sampler.chain)
# lightcurve_param_fit.plot_v_lightcurve_with_fit(SN2018hmx, sampler)


# print(sampler.chain.shape)

# TODO move functions to seperate module files
# TODO write code for extracting tPT via emcee optimization of a0, w0 and p0 (what is m0?)



plt.show()

