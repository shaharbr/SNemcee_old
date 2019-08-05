from lightcurve_lightcurve_fitting import LC
import bolometric
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# TODO compare (plot on top) with and without R band, radius, temp, bolometric lum

# TODO update bolometric table according to new phot data, and also with atlas

lc = LC.read('ascii18hmx.ascii')



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

print(lc)

outpath = '/Users/Shahar/sn2018hmx/SN2018hmx_bolometric_allsources'

t = bolometric.calculate_bolometric(lc, z, outpath, save_table_as='bolometric_tabe_SN2018hmx')

fig = bolometric.plot_bolometric_results(t, save_plot_as='results_fig_allsources.png')

blackbody_data = {'t_from_discovery': np.array(t['MJD']) - 58408.6459, # subtract discovery date
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

blackbody_data.to_csv('blackbody_results.csv')



plt.show()
