import os
os.environ['PYSYN_CDBS'] = os.path.join('..', 'cdbs')
import pandas as pd
import numpy as np
import pysynphot
from matplotlib import pyplot as plt


# data_filepath = os.path.join('results', 'SN2018hmx_lightcurves')
# SN_data = pd.read_csv(data_filepath, usecols=['dmag', 'filter', 'abs_mag', 't_from_discovery'])

# times_list = pd.unique(SN_data['t_from_discovery'])
# +sorted(times_list[times_list > 0.1]* 86400)
# times = sorted((pd.unique(SN_data['t_from_discovery']))* 86400)

times = list(np.arange(0.1, 2.0, 0.1)* 86400)\
        + list(np.arange(2.0, 10.0, 0.5)* 86400)\
        + list(np.arange(10.0, 100.0, 5.0)* 86400)\
        + list(np.arange(100.0, 150.0, 0.5)* 86400)\
        + list(np.arange(150.0, 210.0, 5.0)* 86400)\
        + list(np.arange(210.0, 450.0, 100.0)* 86400)
num_timepoints = len(times)

solar_rad = 6.957 * (10 ** 10)  # cm
z = 0.038 # SN redshift
parsec_dist = z * 4.22 * (10 ** 9) # Parsec

temp_rad_dir = os.path.join('..', 'all_temp_rad_data')
filenames = ['M9.0_Ni0.02_E1.2_Mix1.0_R600_K0.001']
filenames = [str('M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)
                 + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM))
             for Mzams in [9.0]
             for Ni_mass in [0.02]
             for E_final in [1.2]
             for Ni_boundary in [1.0, 3.0, 6.0]
             for R_CSM in [600, 1400, 2200, 3000]
             for K_CSM in [0.001, 50, 100, 150]
             ]


filters = {'u': 'sdss,u', 'g': 'sdss,g', 'r': 'sdss,r', 'i': 'sdss,i', 'z': 'sdss,z',
           'U': 'landolt,u', 'B': 'landolt,b', 'V':'landolt,v', 'R':'landolt,r', 'I':'landolt,i'}
# filters = ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']
# headers = ['Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']

for filename in filenames:
    mag_dict = {'time': times, 'Teff':[], 'PTF_R_AB': np.zeros(num_timepoints)}
    for filtername in ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
        mag_dict[filtername] = []

    temp_table = pd.read_csv(os.path.join(temp_rad_dir, filename, 'T_eff.dat'),
                             names=['time', 'Teff'], sep=r'\s+')
    rad_table = pd.read_csv(os.path.join(temp_rad_dir, filename, 'rad_photo.dat'),
                            names=['time', 'rad'], sep=r'\s+')
    temp_table['time'] = temp_table['time']
    rad_table['time'] = rad_table['time']
    for time in times:
        print(time / 86400)
        temp = np.interp(time, temp_table['time'], temp_table['Teff'])
        mag_dict['Teff'].append(temp)
        print(temp)
        rad = np.interp(time, rad_table['time'], rad_table['rad'])  / solar_rad # in solar radii
        print(rad)
        bb = pysynphot.BlackBody(temp)  # For 1 Solar Radius at 1 kpc
        for filt in filters.keys():
            filter_name = filters[filt]
            bp = pysynphot.ObsBandpass(filter_name)
            obs = pysynphot.Observation(bb, bp)
            if ('sdss' in filter_name) or filter_name == 'o':
                magsystem = 'ABMag'
            else:
                magsystem = 'VegaMag'
            expected_mag = obs.effstim(magsystem) - 2.5*np.log10((rad**2)*((1000.0/parsec_dist)**2))
            mag_dict[filt].append(expected_mag) # Rescaling from the default (1 solar radius at 1000 pc)
    mag_df = pd.DataFrame(mag_dict)
    print(mag_df)
    # plt.plot(mag_df['time']/86400, mag_df['V'] - 36.01, marker='o')
    # plt.gca().invert_yaxis()
    # plt.show()
    # os.mkdir(os.path.join('..', 'all_bb_mag_data', filename))
    mag_df.to_csv(os.path.join('..', 'all_bb_mag_data', filename, 'magnitudes.dat'),
                  sep=' ', header=False, index=False)







