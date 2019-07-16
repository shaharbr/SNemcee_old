import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
import re
import random
import spectra_velocity

# TODO: convert all datetime to MJD

plt.rcParams['font.sans-serif'] = 'Arial'


# initiate SN dictionaries
SN2018hmx = {'name': '2018hmx', 'z': 0.037, 'discovery_date': '2018-10-18 00:00:00', 'spectra': {}, 'expansion_velocities': ''}
SNiPTF14hls = {'name': 'iPTF14hls', 'z': 0.0344, 'discovery_date': '2014-09-22 12:43:00', 'spectra': {}, 'expansion_velocities': ''}

# convert discovery date to datetime format
for SN in [SN2018hmx, SNiPTF14hls]:
    SN['discovery_date'] = pd.to_datetime(SN['discovery_date'])



lines_dict = {'Halpha': {'peak': 6562.81, 'absorption_range': [6200, 6500], 'emission_range': [6400, 6600]},
              'Hbeta': {'peak': 4861, 'absorption_range': [4600, 4800], 'emission_range': [4700, 4900]},
              'FeII 5169': {'peak': 5169, 'absorption_range': [5000, 5100], 'emission_range': [5100, 5500]}
              }

# import expansion velocity file for iPTF14hls
expans_v_iPTF14hls = pd.read_csv('iPTF14hls_expansion_velocity.csv', header=0)
expans_v_iPTF14hls['JD'] = pd.to_datetime(expans_v_iPTF14hls['JD'], unit='D', origin='julian')
expans_v_iPTF14hls.rename(columns={'JD':'datetime', 'Velocity [km/s]':'absorption_mean_velocity', 'Line':'line',
                                   'Velocity_Error [km/s]':'absorption_std_velocity'}, inplace=True)
expans_v_iPTF14hls['name'] = 'iPTF14hls'
expans_v_iPTF14hls['t_from_discovery'] = expans_v_iPTF14hls['datetime'] - SNiPTF14hls['discovery_date']
SNiPTF14hls['expansion_velocities'] = expans_v_iPTF14hls

# impott expansion velocity file for SN1999em
expans_v_sn1999em = pd.read_csv('SN1999em_velocities_leonard_et_al_2002.txt', sep='\s+', header=0)


# import all LCO spectra ascii files for 2018hmx and organize in a dictionary
folder_join = os.path.join
working_dir = r"snexdata_target5025"  # working folder
filenames = os.listdir(working_dir)
# reading and merging
SN2018hmx_spectra = {}
for file in filenames:
    # extract date of spectrum measurment from filename
    date = re.sub('.*hmx|_|-|P60.*|v1|[a-z]|[A-Z]|\..*', '', os.path.basename(file))
    # transform to standard datetime format
    date = pd.to_datetime(date)
    # add as element in dict
    if 'ZTF' in os.path.basename(file):
        SN2018hmx['spectra'][date] = {'df': pd.read_csv(folder_join(working_dir, file), sep=' ', names = ["x", "y", 'dy'], header=180)}
    else:
        SN2018hmx['spectra'][date] = {'df': pd.read_csv(folder_join(working_dir, file), sep=' ', names = ["x", "y"])}
    # SN2018hmx['spectra'][date]['df']['x'] = SN2018hmx['spectra'][date]['df']['x'].astype('float')
    # SN2018hmx['spectra'][date]['df']['y'] = SN2018hmx['spectra'][date]['df']['y'].astype('float')
    if 'redblu' in os.path.basename(file):
        SN2018hmx['spectra'][date]['telescope'] = 'LCO'
    elif 'ZTF' in os.path.basename(file):
        SN2018hmx['spectra'][date]['telescope'] = 'ZTF'
    elif 'HET' in os.path.basename(file):
        SN2018hmx['spectra'][date]['telescope'] = 'HET'
    else:
        SN2018hmx['spectra'][date]['telescope'] = 'ND'

# import ZTF spectra from Weizmann transient name server as ascii file and add it to the dictionary
# filename = 'tns_2018hmx_2018-11-06_09-27-00_P60_SED-Machine_ZTF.ascii'
# date = re.sub('tns.*hmx_|_P.*', '', filename)
# date = re.sub('_', ' ', date)
# transform to standard datetime format
# date = pd.to_datetime(date)
# add as element in dict
# SN2018hmx['spectra'][date] = {'df': pd.read_csv(filename, sep=' ', names = ["x", "y", 'dy'], header=180)}
# correct wavelengths for redshift

SN2018hmx = spectra_velocity.add_time_from_discovery(SN2018hmx)


SN2018hmx = spectra_velocity.correct_redshift(SN2018hmx)


SN2018hmx = spectra_velocity.normalize_spectra(SN2018hmx)

spectra_velocity.plot_overlay_spectra(SN2018hmx, lines_dict)

SN2018hmx['spectra'] = spectra_velocity.fit_Pcygni_curves(SN2018hmx['spectra'], lines_dict, fixed_curve_range=False, number_curves=3)

spectra_velocity.plot_stacked_spectra(SN2018hmx, lines_dict)

spectra_velocity.plot_stacked_spectra(SN2018hmx, lines_dict, plot_curve_fits=True)

SN2018hmx['spectra'] = spectra_velocity.add_expansion_velocity(SN2018hmx['spectra'], lines_dict)


SN2018hmx['expansion_velocities'] = spectra_velocity.make_velocity_df(SN2018hmx, lines_dict)

expans_v_sn1999em = expans_v_sn1999em.loc[expans_v_sn1999em['day'] < 96]


spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SNiPTF14hls['expansion_velocities']], 'absorption')
# plot_expansion_velocities([SN2018hmx['expansion_velocities']], 'emission')

# SN2018hmx['expansion_velocities'].to_csv('sN2018hmx_expansion_velocities.csv')

plt.show()

