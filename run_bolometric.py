from lightcurve_fitting import lightcurve, models, fitting, bolometric
from pkg_resources import resource_filename
from IPython.display import display, Math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os




# TODO compare (plot on top) with and without R band, radius, temp, bolometric lum

# TODO update bolometric table according to new phot data, and also with atlas

SN = '04et'


filename = os.path.join('data', 'ascii'+SN+'.ascii')
correction_data = pd.read_csv(os.path.join('data', 'SN20'+SN+'_correction.csv'))
discovery_date = float(correction_data['discovery_date'])

extinction = {}
filters = ['U', 'B', 'V', 'R', 'I', 'u', 'g', 'r', 'i', 'z', 'J', 'H', 'K', 'L', 'G', 'o']
for filter in filters:
    if filter == 'z':
        extinction[filter] = correction_data['z_sdss'][0]
    else:
        extinction[filter] = correction_data[filter][0]


lc = lightcurve.LC.read(filename)
lc.meta['dm'] = float(correction_data['dm'])
lc.meta['extinction'] = extinction
z = float(correction_data['z'])
print(lc)
outpath = SN+'_BVgri_bolometric_allsources'


t = bolometric.calculate_bolometric(lc, z, outpath, colors=['B-V', 'g-r', 'r-i'], save_table_as='bolometric_table_'+SN+'_BVgri')
print(t)
fig1 = bolometric.plot_bolometric_results(t, save_plot_as=SN+'_bolometric.png')
fig2 = bolometric.plot_color_curves(t)

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

blackbody_data.to_csv(os.path.join('results', 'blackbody_results_SN20'+SN+'_BVgri.csv'))


plt.show()
