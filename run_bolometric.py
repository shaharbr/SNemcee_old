from lightcurve_lightcurve_fitting import LC
import bolometric
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# TODO compare (plot on top) with and without R band, radius, temp, bolometric lum

# TODO update bolometric table according to new phot data, and also with atlas

filename = os.path.join('data', 'ascii18hmx.ascii')

SN = '18aad'

if SN == '18hmx':
    filename = os.path.join('data', 'ascii18hmx.ascii')
    discovery_date = 58408.646
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
    z = 0.038
    dm = 36.01

if SN == '18aad':
    filename = os.path.join('data', 'ascii18aad.ascii')
    discovery_date = 58182.7996
    extinction = {
        'U': 0,
        'B': 0,
        'V': 0,
        'R': 0,
        'I': 0,
        'u': 0,
        'g': 0,
        'r': 0,
        'i': 0,
        'z': 0,
        'J': 0,
        'H': 0,
        'K': 0,
        'L': 0,
        'G': 0,
        'o': 0
    }
    z = 0.024
    dm = 35.178


lc = LC.read(filename)
lc.calcAbsMag(dm=dm, extinction=extinction)
lc.calcLum()

print(lc)


outpath = SN+'_BVgri_bolometric_allsources'
table_path = os.path.join(outpath, 'bolometric_tabe_'+SN+'BVgri')
fig_path = os.path.join(outpath, 'results_fig_allsources.png')

# note: when astropy is version >=4, there is a problem with trying to convert a dimentionless value here
# (works with 3.2.3)
t = bolometric.calculate_bolometric(lc, z, outpath, save_table_as='bolometric_tabe_'+SN+'_BVgri')
fig = bolometric.plot_bolometric_results(t, save_plot_as='results_fig_allsources.png')

blackbody_data = {'t_from_discovery': np.array(t['MJD']) - discovery_date, # subtract discovery date
                  'temp': np.array(t['temp_mcmc']),
                  'dtemp0': np.array(t['dtemp0']),
                  'dtemp1': np.array(t['dtemp1']),
                  'radius': np.array(t['radius_mcmc']) * 695510000, # convert "1000 solar radius" units to km
                  'dradius0': np.array(t['dradius0']) * 695510000,
                  'dradius1': np.array(t['dradius1']) * 695510000,
                  'Lum': np.array(t['L_mcmc']),
                  'dLum0': np.array(t['dL_mcmc0']),
                  'dLum1': np.array(t['dL_mcmc1']),
                  }

blackbody_data = pd.DataFrame.from_dict(blackbody_data)

blackbody_data.to_csv(os.path.join('results', 'blackbody_results_'+SN+'_BVgri.csv'))


plt.show()
