from lightcurve_lightcurve_fitting import LC
import bolometric
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# TODO compare (plot on top) with and without R band, radius, temp, bolometric lum

# TODO update bolometric table according to new phot data, and also with atlas

lc = LC.read(r'data\ascii18aad_BVgri.ascii')

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
lc.calcAbsMag(dm=35.21, extinction=extinction)
lc.calcLum()

print(lc)

outpath = '/Users/Shahar/sn2018hmx/SN2018aad_BVgri_bolometric_allsources'

t = bolometric.calculate_bolometric(lc, z, outpath, save_table_as='bolometric_tabe_SN2018aad')

# fig = bolometric.plot_bolometric_results(t, save_plot_as='results_fig_allsources.png')
fig = bolometric.plot_bolometric_results(t)

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

blackbody_data.to_csv(r'results\blackbody_results_18aad_BVgri.csv')



plt.show()
