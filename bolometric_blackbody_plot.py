from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import Ni_mass



'''
Fitting the Planck function using an MCMC routine. This is slower, depending on how many walkers (nwalkers) and steps (burnin_steps and steps) you use, but gives more robust uncertainties. The columns temp_mcmc, radius_mcmc, dtemp0, dtemp1, dradius0, and dradius1 come from this fit. My convention for non-Gaussian uncertainties is that 0 is the lower uncertainty and 1 is the upper uncertainty.
Integrating the Planck function between $U$ and $I$ band (observed) gives L_mcmc, dL_mcmc0, and dL_mcmc1.
'''

SN1987A_bolometric = pd.read_csv(r'data\bersten_1987a_bol.csv', names=['t_from_discovery', 'Lum'])
SN1987A_bolometric['Lum'] = 10 ** SN1987A_bolometric['Lum']


# SN87A_line, SN87A_sigma = Ni_mass.calc_87A_slope(SN1987A_bolometric)

# fit
# SN87A_tail = SN1987A_bolometric.loc[SN1987A_bolometric['t_from_discovery'] > 132]
# x_87A = SN87A_tail['t_from_discovery']
# y_87A = np.log10(SN87A_tail['Lum'])
# SN87A_line = np.polyfit(x_87A, y_87A, deg=1)
# SN87A_line, SN87A_cov = np.polyfit(x_87A, y_87A, deg=1, cov=True)
# SN87A_sigma = np.sqrt(np.diag(SN87A_cov))
# SN87A_sigma = SN87A_sigma[0]



# TODO fiish turning everythong to For loop

# TODO data should load into dicts and then the loop runs on the dict by name of SN
blackbody_data = {'SN2018hmx': [], 'SN2018aad': []}
velocity_data = {'SN2018hmx': [], 'SN2018aad': []}
# blackbody_data['SN2018hmx'] = blackbody_data['SN2018hmx'].loc[blackbody_data['SN2018hmx'] > ]
blackbody_data['SN2018hmx'] = pd.read_csv(r'results\blackbody_results_18hmx_BVgi.csv')

velocity_data['SN2018hmx'] = pd.read_csv(r'results\sN2018hmx_expansion_velocities.csv')
print(blackbody_data['SN2018hmx'])
# exit()
blackbody_data['SN2018aad'] = pd.read_csv(r'results\blackbody_results_18aad_BVgri.csv')
velocity_data['SN2018aad'] = pd.read_csv(r'results\SN2018aad_expansion_velocities.csv')

# TODO something wrong here with the Ni mass calculation of 18hmx

# convert watt to erg/s
def convert_to_erg_s(blackbody_data):
    blackbody_data['Lum'] = blackbody_data['Lum'] * 10 ** 7
    blackbody_data['dLum0'] = blackbody_data['dLum0'] * 10 ** 7
    blackbody_data['dLum1'] = blackbody_data['dLum1'] * 10 ** 7
    return blackbody_data



for SN in ['SN2018hmx']:
    blackbody_data[SN] = convert_to_erg_s(blackbody_data[SN])


for SN in ['SN2018hmx']:
    print(SN)
    velocity_data[SN] = velocity_data[SN].loc[velocity_data[SN]['line'] == 'FeII 5169'].reset_index()
    velocity_data[SN]['vt'] = velocity_data[SN]['t_from_discovery'] * velocity_data[SN]['absorption_mean_velocity'] * 86400 #vt and multiplied by seconds in a day (v in km/s)
    velocity_data[SN]['d_vt'] = velocity_data[SN]['t_from_discovery'] * velocity_data[SN]['absorption_std_velocity'] * 86400

    plt.figure()
    plt.errorbar(x=blackbody_data[SN]['t_from_discovery'],
                y=blackbody_data[SN]['radius'],
                yerr=[blackbody_data[SN]['dradius0'],
                      blackbody_data[SN]['dradius1']],
                label=SN + ' blackbody radius',
                marker='o',
                fillstyle='full',
                linestyle='None',
                color='green')

    plt.errorbar(x=velocity_data[SN]['t_from_discovery'],
                y=velocity_data[SN]['vt'],
                yerr=velocity_data[SN]['d_vt'],
                label=SN + ' Fe velocity * t',
                marker='o',
                fillstyle='full',
                linestyle='None',
                color='blue')


    # plt.yscale('log')
    plt.ylabel('radius (km)')
    plt.xlabel('time after discovery (days)')
    plt.title(SN+' - blackbody radius vs Fe v*t\n', fontsize=14)
    plt.ylim()
    # ax1.legend(ncol=2)
    plt.legend()



    plt.figure()
    plt.errorbar(x=blackbody_data[SN]['t_from_discovery'],
                y=blackbody_data[SN]['Lum'],
                yerr=[blackbody_data[SN]['dLum0'],
                      blackbody_data[SN]['dLum1']],
                label=SN,
                marker='o',
                fillstyle='full',
                linestyle='None',
                color='black')

    plt.errorbar(x=SN1987A_bolometric['t_from_discovery'],
                y=SN1987A_bolometric['Lum'],
                label='SN 1987A',
                marker='o',
                fillstyle='full',
                linestyle='None',
                color='gray')

    plt.yscale('log')
    plt.xlim(left=0)
    plt.ylabel('bolometric luminosity (erg/s)')
    plt.xlabel('time after discovery (days)')
    plt.title(SN + ' - bolometric luminosity \n', fontsize=14)
    plt.legend()

    # fit






    # y_fit_87A = np.power(10, SN87A_line[1] + x * SN87A_line[0])


    x_fit_SN, y_fit_SN, y_fit_SN_sigma = Ni_mass.fit_to_log_slope(blackbody_data[SN])
    x_fit_87A, y_fit_87A, y_fit_87A_sigma = Ni_mass.fit_to_log_slope(SN1987A_bolometric)
    # plt.figure()
    # plt.plot(x_fit_SN, y_fit_SN_sigma)
    # plt.figure()
    # plt.plot(x_fit_87A,y_fit_87A_sigma)


    plt.errorbar(x=x_fit_SN, y=y_fit_SN, yerr=y_fit_SN_sigma)
    plt.errorbar(x=x_fit_87A, y=y_fit_87A, yerr=y_fit_87A_sigma)




    sampler = Ni_mass.SN_lightcurve_params(blackbody_data[SN])
    Ni_mass.chain_plots(sampler.chain)
    Ni_vec = Ni_mass.plot_v_lightcurve_with_fit(blackbody_data[SN], sampler)

    Ni_df = pd.DataFrame({'t_from_disovery': blackbody_data[SN]['t_from_discovery'], 'Ni_mass': Ni_vec})
    Ni_df.to_csv(r'results\Ni_mass_'+SN+'_BVgri.csv')

    Ni_results = Ni_mass.get_Ni_results_dict(sampler)

    Ni_results['Ni_87A'], Ni_results['Ni_87A_sigma'] =\
        Ni_mass.Ni_by_87A_slope(y_fit_SN, y_fit_SN_sigma, y_fit_87A, y_fit_87A_sigma)
    print(Ni_results)
    param_results = pd.DataFrame(Ni_results, index=[0])
    param_results.to_csv(r'results\Ni_results_'+SN+'_BVgri.csv')





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


# blackbody_radius_zero(blackbody_data)






plt.show()
