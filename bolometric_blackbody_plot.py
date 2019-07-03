from lightcurve import LC
import bolometric
from matplotlib import pyplot as plt
import pandas as pd
import data_import
import numpy as np



'''
Fitting the Planck function using an MCMC routine. This is slower, depending on how many walkers (nwalkers) and steps (burnin_steps and steps) you use, but gives more robust uncertainties. The columns temp_mcmc, radius_mcmc, dtemp0, dtemp1, dradius0, and dradius1 come from this fit. My convention for non-Gaussian uncertainties is that 0 is the lower uncertainty and 1 is the upper uncertainty.
Integrating the Planck function between $U$ and $I$ band (observed) gives L_mcmc, dL_mcmc0, and dL_mcmc1.
'''


expansion_v = pd.read_csv('../sn2018hmx/sN2018hmx_expansion_velocities.csv')
expansion_v = expansion_v.loc[expansion_v['line'] == 'FeII 5169'].reset_index()
expansion_v['datetime'] = data_import.convert_to_mjd(expansion_v['datetime'], from_datetime=True)
expansion_v['t_from_discovery'] = expansion_v['datetime'] - 58408.6459
expansion_v['vt'] = expansion_v['t_from_discovery'] * expansion_v['absorption_mean_velocity'] * 86400 #vt and multiplied by seconds in a day (v in km/s)
expansion_v['d_vt'] = expansion_v['t_from_discovery'] * expansion_v['absorption_std_velocity'] * 86400



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
