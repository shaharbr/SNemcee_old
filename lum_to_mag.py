import os
os.environ['PYSYN_CDBS'] = os.path.join('..', 'cdbs')
import pandas as pd
import numpy as np
import pysynphot
from matplotlib import pyplot as plt



dist_mod = 36.01
z = 0.038 # SN redshift
parsec_dist = z * 4.22 * (10 ** 9) # Parsec
solar_rad = 6.957 * (10 ** 10)  # cm

colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}

times = list(np.arange(0.1, 2.0, 0.1)* 86400)\
        + list(np.arange(2.0, 10.0, 0.5)* 86400)\
        + list(np.arange(10.0, 100.0, 5.0)* 86400)\
        + list(np.arange(100.0, 150.0, 0.5)* 86400)\
        + list(np.arange(150.0, 210.0, 5.0)* 86400)\
        + list(np.arange(210.0, 450.0, 100.0)* 86400)
# num_timepoints = len(times)


temp_rad_dir = os.path.join('..', 'all_temp_rad_data')
mag_dir = os.path.join('..', 'all_mag_data')
lum_dir = os.path.join('..', 'all_lum_data')

filenames = [str('M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)
                 + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM))
             for Mzams in [15.0]
             for Ni_mass in [0.07]
             for E_final in [1.2]
             for Ni_boundary in [3.0]
             for R_CSM in [600]
             for K_CSM in [0.001]
             ]


filters = {'u': 'sdss,u', 'g': 'sdss,g', 'r': 'sdss,r', 'i': 'sdss,i', 'z': 'sdss,z',
           'U': 'landolt,u', 'B': 'landolt,b', 'V':'landolt,v', 'R':'landolt,r', 'I':'landolt,i'}




# filters = ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']
# headers = ['Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']


for filename in filenames:
    temp_table = pd.read_csv(os.path.join(temp_rad_dir, filename, 'T_eff.dat'),
                             names=['time', 'Teff'], sep=r'\s+')
    rad_table = pd.read_csv(os.path.join(temp_rad_dir, filename, 'rad_photo.dat'),
                            names=['time', 'rad'], sep=r'\s+')
    lum_table = pd.read_csv(os.path.join(lum_dir, filename, 'lum_observed.dat'),
                            names=['time', 'lum'], sep=r'\s+')

    plt.figure()
    plt.plot(temp_table['time'] / 86400, temp_table['Teff'])
    plt.title('Effective temperature over time\n'
              + filename)
    plt.xlabel('time (days)')
    plt.ylabel('Teff (K)')
    plt.ylim(0, 10000)
    plt.tight_layout()
    plt.savefig(filename + 'temp_zoom.png')

    plt.figure()
    plt.plot(rad_table['time'] / 86400, rad_table['rad'])
    plt.title('Radius  over time\n'
              + filename)
    plt.xlabel('time (days)')
    plt.ylabel('radius (R_sol)')
    plt.tight_layout()
    plt.savefig(filename+'radius.png')

    plt.figure()
    plt.plot(lum_table['time'] / 86400, lum_table['lum'])
    plt.title('Bolometric luminosity  over time\n'
              + filename)
    plt.xlabel('time (days)')
    plt.ylabel('Luminosity (erg/s)')
    plt.ylim(1e+40, 0.4e+43)
    plt.tight_layout()
    plt.yscale('log')
    plt.savefig(filename + 'lum.png')


snecmag = pd.read_csv(os.path.join(mag_dir, filename, 'magnitudes.dat'),
                           names=['time', 'Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'],
                           sep=r'\s+')
snecmag = snecmag.abs()


for filename in filenames:
    temp_table = pd.read_csv(os.path.join(temp_rad_dir, filename, 'T_eff.dat'),
                             names=['time', 'Teff'], sep=r'\s+')
    rad_table = pd.read_csv(os.path.join(temp_rad_dir, filename, 'rad_photo.dat'),
                            names=['time', 'rad'], sep=r'\s+')
    temp_table['time'] = temp_table['time']
    rad_table['time'] = rad_table['time']

    times = list(temp_table['time'])
    times_smaller = times[1:100:]+times[100::10]
    num_timepoints = len(times_smaller)
    print(num_timepoints)

    mag_dict = {'time': times_smaller, 'Teff': [], 'PTF_R_AB': np.zeros(num_timepoints)}
    for filtername in ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
        mag_dict[filtername] = []

    for time in times_smaller:
        print(time)
        # temp = np.interp(time, temp_table['time'], temp_table['Teff'])
        temp = float(temp_table.loc[temp_table['time'] == time]['Teff'])
        mag_dict['Teff'].append(temp)
        print(temp)
        # rad = np.interp(time, rad_table['time'], rad_table['rad'])  / solar_rad # in solar radii
        rad = float(rad_table.loc[rad_table['time'] == time]['rad'] / solar_rad)
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
            expected_mag = -(expected_mag - dist_mod)# remove distance modulus and get positive value
            mag_dict[filt].append(expected_mag) # Rescaling from the default (1 solar radius at 1000 pc)
    mag_df = pd.DataFrame(mag_dict)
    print(mag_df)
    # plt.plot(mag_df['time']/86400, mag_df['V'] - 36.01, marker='o')
    # plt.gca().invert_yaxis()
    # plt.show()
    # os.mkdir(os.path.join('..', 'all_bb_mag_data', filename))
    # mag_df.to_csv(os.path.join('..', 'all_bb_mag_data', filename, 'magnitudes.dat'),
    #               sep=' ', header=False, index=False)


    titles =['Pysynphot', 'SNEC']
    i = 0
    for magtable in [mag_df, snecmag]:
        f_fit, ax = plt.subplots()
        for filt in ['g', 'r', 'i', 'U', 'B', 'V']:
            ax.plot(magtable['time'] / 86400, -magtable[filt], color=colors[filt], label=filt)
        ax.set_title(titles[i]+' Magnitudes over time\n'
                  + filename)
        ax.set_xlabel('time (days)')
        ax.set_ylabel('Absolute magnitude')
        ax.set_ylim(-20, 30)
        ax.legend()
        plt.gca().invert_yaxis()
        plt.tight_layout()
        f_fit.savefig(filename + titles[i]+'_mag.png')
        i +=1
    i = 0
    for magtable in [mag_df, snecmag]:
        f_fit, ax = plt.subplots()
        for filt in ['g', 'r', 'i', 'U', 'B', 'V']:
            ax.plot(magtable['time'] / 86400, -magtable[filt], color=colors[filt], label=filt)
        ax.set_title(titles[i]+' Magnitudes over time\n'
                  + filename)
        ax.set_xlabel('time (days)')
        ax.set_ylabel('Absolute magnitude')
        ax.set_ylim(-18, -11)
        ax.legend()
        plt.gca().invert_yaxis()
        plt.tight_layout()
        f_fit.savefig(filename + titles[i]+'_mag_zoom.png')
        i +=1



plt.show()




