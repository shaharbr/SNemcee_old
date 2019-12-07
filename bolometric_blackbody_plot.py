from matplotlib import pyplot as plt
import pandas as pd
import data_import_lightcurve_fitting as data_import
import numpy as np
import Ni_mass


'''
Fitting the Planck function using an MCMC routine. This is slower, depending on how many walkers (nwalkers) and steps (burnin_steps and steps) you use, but gives more robust uncertainties. The columns temp_mcmc, radius_mcmc, dtemp0, dtemp1, dradius0, and dradius1 come from this fit. My convention for non-Gaussian uncertainties is that 0 is the lower uncertainty and 1 is the upper uncertainty.
Integrating the Planck function between $U$ and $I$ band (observed) gives L_mcmc, dL_mcmc0, and dL_mcmc1.
'''

# SN14hls_radius_temp = pd.read_csv(r'results\14hls_blackbody_radius_temp.csv')
# SN14hls_bolometric = pd.read_csv(r'results\14hls_bolometric_luminosity.csv')
SN1987A_bolometric = pd.read_csv(r'data\bersten_1987a_bol.csv', names=['t_from_discovery', 'Lum'])
SN1987A_bolometric['Lum'] = 10 ** SN1987A_bolometric['Lum']




# get d uncertainties
# SN14hls_radius_temp['Radius_Low [cm]'] = SN14hls_radius_temp['Radius [cm]'] - SN14hls_radius_temp['Radius_Low [cm]']
# SN14hls_radius_temp['Radius_High [cm]'] = SN14hls_radius_temp['Radius_High [cm]'] - SN14hls_radius_temp['Radius [cm]']

# SN14hls_bolometric['Luminosity_Low [erg/s]'] = SN14hls_bolometric['Luminosity [erg/s]'] - SN14hls_bolometric['Luminosity_Low [erg/s]']
# SN14hls_bolometric['Luminosity_High [erg/s]'] = SN14hls_bolometric['Luminosity_High [erg/s]'] - SN14hls_bolometric['Luminosity [erg/s]']




# convert JD to MJD
# for df in [SN14hls_radius_temp, SN14hls_bolometric]:
#     df['JD'] = data_import.convert_to_mjd(df['JD'])
#     df.rename(columns={'JD': 'mjd'}, inplace=True)
#     df['t_from_discovery'] = df['mjd'] - 56922.52986111




# convert cm to km
# SN14hls_radius_temp['Radius [cm]'] = SN14hls_radius_temp['Radius [cm]'] / 10**5
# SN14hls_radius_temp['Radius_Low [cm]'] = SN14hls_radius_temp['Radius_Low [cm]'] / 10**5
# SN14hls_radius_temp['Radius_High [cm]'] = SN14hls_radius_temp['Radius_High [cm]'] / 100**5
# SN14hls_radius_temp.rename(columns={'Radius [cm]': 'Radius [km]', 'Radius_Low [cm]': 'Radius_Low [km]', 'Radius_High [cm]': 'Radius_High [km]'}, inplace=True)




# SN14hls_vt = SN14hls_radius_temp.loc[SN14hls_radius_temp['Radius_Type'] == 'Fe II v*t']
# SN14hls_blackbody = SN14hls_radius_temp.loc[SN14hls_radius_temp['Radius_Type'] == 'Blackbody']


blackbody_data = pd.read_csv(r'results\blackbody_results_18aad_BVgri.csv')


# convert watt to erg/s
blackbody_data['Lum'] = blackbody_data['Lum'] * 10**7
blackbody_data['dLum0'] = blackbody_data['dLum0'] * 10**7
blackbody_data['dLum1'] = blackbody_data['dLum1'] * 10**7






# expansion_v = pd.read_csv(r'results\sN2018hmx_expansion_velocities.csv')
# expansion_v = expansion_v.loc[expansion_v['line'] == 'FeII 5169'].reset_index()
# expansion_v['vt'] = expansion_v['t_from_discovery'] * expansion_v['absorption_mean_velocity'] * 86400 #vt and multiplied by seconds in a day (v in km/s)
# expansion_v['d_vt'] = expansion_v['t_from_discovery'] * expansion_v['absorption_std_velocity'] * 86400



fig1, ax1 = plt.subplots()
ax1.errorbar(x=blackbody_data['t_from_discovery'],
            y=blackbody_data['radius'],
            yerr=[blackbody_data['dradius0'],
                  blackbody_data['dradius1']],
            label='18aad blackbody radius',
            marker='o',
            fillstyle='full',
            linestyle='None',
            color='green')

# ax1.errorbar(x=expansion_v['t_from_discovery'],
#             y=expansion_v['vt'],
#             yerr=expansion_v['d_vt'],
#             label='18hmx Fe velocity * t',
#             marker='o',
#             fillstyle='full',
#             linestyle='None',
#             color='blue')

# ax1.errorbar(x=SN14hls_vt['t_from_discovery'],
#             y=SN14hls_vt['Radius [km]'],
#             yerr=[SN14hls_vt['Radius_Low [km]'],
#                   SN14hls_vt['Radius_High [km]']],
#             label='14hls blackbody radius',
#             marker='o',
#             fillstyle='none',
#             linestyle='None',
#             color='lightblue')

# ax1.errorbar(x=SN14hls_blackbody['t_from_discovery'],
#             y=SN14hls_blackbody['Radius [km]'],
#             yerr=[SN14hls_blackbody['Radius_Low [km]'],
#                   SN14hls_blackbody['Radius_High [km]']],
#             label='14hls Fe velocity * t',
#             marker='o',
#             fillstyle='none',
#             linestyle='None',
#             color='lightgreen')
# plt.yscale('log')
ax1.set_ylabel('radius (km)')
ax1.set_xlabel('time after discovery (days)')
ax1.set_title('SN2018aad - blackbody radius vs Fe v*t\n', fontsize=14)
ax1.set_ylim()
# ax1.legend(ncol=2)
ax1.legend()



fig2, ax2 = plt.subplots()
ax2.errorbar(x=blackbody_data['t_from_discovery'],
            y=blackbody_data['Lum'],
            yerr=[blackbody_data['dLum0'],
                  blackbody_data['dLum1']],
            label='SN 2018aad',
            marker='o',
            fillstyle='full',
            linestyle='None',
            color='black')

ax2.errorbar(x=SN1987A_bolometric['t_from_discovery'],
            y=SN1987A_bolometric['Lum'],
            label='SN 1987A',
            marker='o',
            fillstyle='full',
            linestyle='None',
            color='gray')




# ax2.errorbar(x=SN14hls_bolometric['t_from_discovery'],
#             y=SN14hls_bolometric['Luminosity [erg/s]']
#             yerr=[SN14hls_bolometric['Luminosity_Low [erg/s]']
#                   SN14hls_bolometric['Luminosity_High [erg/s]']],
#             label='14hls',
#             marker='o',
#             fillstyle='none',
#             linestyle='None',
#             color='lightgreen')


ax2.set_yscale('log')
ax2.set_xlim(left=0)
ax2.set_ylabel('bolometric luminosity (erg/s)')
ax2.set_xlabel('time after discovery (days)')
ax2.set_title('SN2018aad - bolometric luminosity \n', fontsize=14)
ax2.legend()




# fit
SN87A_tail = SN1987A_bolometric.loc[SN1987A_bolometric['t_from_discovery'] > 132]
x = SN87A_tail['t_from_discovery']
y = np.log10(SN87A_tail['Lum'])
SN87A_line = np.polyfit(x, y, deg=1)
SN87A_line, SN87A_cov = np.polyfit(x, y, deg=1, cov=True)
SN87A_sigma = np.sqrt(np.diag(SN87A_cov))
SN87A_sigma = SN87A_sigma[0]
y_fit = np.power(10, SN87A_line[1] + x * SN87A_line[0])
ax2.errorbar(x=x, y=y_fit, yerr=SN87A_sigma)

# fit
SN18aad_tail = blackbody_data.loc[blackbody_data['t_from_discovery'] > 132]
x = SN18aad_tail['t_from_discovery']
y = np.log10(SN18aad_tail['Lum'])
SN18aad_line = np.polyfit(x, y, deg=1)
SN18aad_line, SN18aad_cov = np.polyfit(x, y, deg=1, cov=True)
SN18aad_sigma = np.sqrt(np.diag(SN18aad_cov))
SN18aad_sigma = SN18aad_sigma[0]
y_fit = np.power(10, SN18aad_line[1] + x * SN18aad_line[0])
ax2.errorbar(x=x, y=y_fit, yerr=SN18aad_sigma)


def log_slope_ratio(slope1, slope2):
    return np.power(10, slope1 - slope2)

day = 135
print('day 135', 0.075 * (np.power(10, SN18aad_line[1] + day * SN18aad_line[0]) / np.power(10, SN87A_line[1] + day * SN87A_line[0])))

day = 145
print('day 145', 0.075 * (np.power(10, SN18aad_line[1] + day * SN18aad_line[0]) / np.power(10, SN87A_line[1] + day * SN87A_line[0])))

day = 155
print('day 155', 0.075 * (np.power(10, SN18aad_line[1] + day * SN18aad_line[0]) / np.power(10, SN87A_line[1] + day * SN87A_line[0])))

day = 165
print('day 165', 0.075 * (np.power(10, SN18aad_line[1] + day * SN18aad_line[0]) / np.power(10, SN87A_line[1] + day * SN87A_line[0])))

day = 200
print('day 200', 0.075 * (np.power(10, SN18aad_line[1] + day * SN18aad_line[0]) / np.power(10, SN87A_line[1] + day * SN87A_line[0])))




sampler = Ni_mass.SN_lightcurve_params(blackbody_data)
Ni_mass.chain_plots(sampler.chain)
Ni_vec = Ni_mass.plot_v_lightcurve_with_fit(blackbody_data, sampler)

Ni_df = pd.DataFrame({'t_from_disovery': blackbody_data['t_from_discovery'], 'Ni_mass': Ni_vec})
Ni_df.to_csv(r'results\Ni_mass_SN2018aad_BVgri.csv')

Ni_results = Ni_mass.get_Ni_results_dict(sampler)
print(Ni_results)
param_results = pd.DataFrame(Ni_results, index=[0])
param_results.to_csv(r'results\SN2018aad_Ni_results_BVgri.csv')



print(sampler.chain.shape)


def blackbody_radius_zero(SN_df):
    x = SN_df['t_from_discovery'][0:7]
    y = SN_df['radius'][0:7]
    plt.figure()
    plt.plot(x, y, marker='o')
    f = np.poly1d(np.polyfit(x, y, deg=1))
    fit_x = range(100)
    plt.plot(fit_x, f[1] * fit_x + f[0])
    plt.axhline(y=0, color='k', linestyle='--')
    radius_0_t = - f[0] / f[1]
    plt.xlabel('time from discovery (days)')
    plt.ylabel('blackbody radius (km)')
    print(radius_0_t)


blackbody_radius_zero(blackbody_data)






plt.show()
