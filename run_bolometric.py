from lightcurve import LC
import bolometric
from matplotlib import pyplot as plt
import pandas as pd
import data_import
import numpy as np


expansion_v = pd.read_csv('../sn2018hmx/sN2018hmx_expansion_velocities.csv')
expansion_v = expansion_v.loc[expansion_v['line'] == 'FeII 5169'].reset_index()
expansion_v['datetime'] = data_import.convert_to_mjd(expansion_v['datetime'], from_datetime=True)
expansion_v['t_from_discovery'] = expansion_v['datetime'] - 58408.6459
expansion_v['vt'] = expansion_v['t_from_discovery'] * expansion_v['absorption_mean_velocity'] * 86400 #vt and multiplied by seconds in a day (v in km/s)
expansion_v['d_vt'] = expansion_v['t_from_discovery'] * expansion_v['absorption_std_velocity'] * 86400
print(expansion_v['t_from_discovery'])
print(expansion_v['absorption_mean_velocity'])
print(expansion_v['vt'])
print(expansion_v['d_vt'])




lc = LC.read('example/SN2018hmx.table')
# lc = LC.read('example/SN2016bkv.table')

'''
extinction = {
 'U': 0.069,
 'B': 0.061,
 'g': 0.055,
 'V': 0.045,
 '0': 0.035,
 'r': 0.038,
 'R': 0.035,
 'i': 0.028,
 'I': 0.025,
}
z = 0.002
lc.calcAbsMag(dm=30.79, extinction=extinction)
lc.calcLum()

'''


extinction = {
    'U': 0.206,
    'B': 0.173,
    'V': 0.131,
    'R': 0.103,
    'I': 0.072,
    'u': 0.202,
    'g': 0.157,
    'r': 0.109,
    'i': 0.081,
    'z': 0.060,
    'J': 0.034,
    'H': 0.021,
    'K': 0.014,
    'L': 0.007,
    'G': 0,
    'o': 0
}

z = 0.037
lc.calcAbsMag(dm=36.06, extinction=extinction)
lc.calcLum()


outpath = '/Users/Shahar/sn2018hmx/SN2018hmx_bolometric_allsources'
# outpath = '/Users/Shahar/sn2018hmx/SN2016bkv_bolometric'

t = bolometric.calculate_bolometric(lc, z, outpath, save_table_as='bolometric_tabe_SN2018hmx')

fig = bolometric.plot_bolometric_results(t, save_plot_as='results_fig_allsources.png')

print(t)

print(list(t['MJD']))
print(list(t['temp']))
print(list(t['radius']))

print(t['temp'] * t['radius'])
print(list(t['temp'] * t['radius']))

blackbody_data = {'t_from_discovery': np.array(t['MJD']) - 58408.6459, # subtract discovery date
                  # 'temp': np.array(t['temp']),
                  # 'dtemp': np.array(t['dtemp']),
                  'radius': np.array(t['radius']) * 695510000, # convert "1000 solar radius" units to km
                  'dradius': np.array(t['dradius']) * 695510000
                  }





fig1, ax1 = plt.subplots()
ax1.errorbar(x=blackbody_data['t_from_discovery'],
            y=blackbody_data['radius'] / 10000000000,
            yerr=blackbody_data['dradius'] / 10000000000,
            label='blackbody radius',
            marker='o',
            linestyle='None',
            color='blue')
ax1.errorbar(x=expansion_v['t_from_discovery'],
            y=expansion_v['vt'] / 10000000000,
            yerr=expansion_v['d_vt'] / 10000000000,
            label='Fe velocity * t',
            marker='o',
            linestyle='None',
            color='red')
ax1.set_ylabel('radius (10^10 km)')
ax1.set_xlabel('time after discovery (days)')
ax1.set_title('SN2018 - blackbody radius vs Fe velocity * time')
ax1.legend()

plt.show()

