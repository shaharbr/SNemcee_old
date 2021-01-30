import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import mcmc_snec
import re
import snec_result_interpolator_lum as interp_lum
import snec_result_interpolator_veloc as interp_veloc



result_t = '2021-01-18_00-11-37_lum_twostep_flatpriors_SN2017eaw'
result_dir = os.path.join('mcmc_results', result_t)

flat_sampler_path = os.path.join(result_dir, 'flat_sampler.csv')
run_params_path = os.path.join(result_dir, 'run_parameters.csv')
output_path = os.path.join(result_dir, 'run_parameters.csv')

run_params = pd.read_csv(run_params_path, index_col=0).T


SN_name = run_params.iloc[0]['SN_name']
n_walkers = int(run_params.iloc[0]['n_walkers'])
n_steps = int(run_params.iloc[0]['n_steps'])


distances = pd.read_csv(os.path.join('results', 'distances.csv'))
distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
sigma_S = 2 * distance_err / (distance)


Mzams_range = [9.0, 11.0, 13.0, 15.0, 17.0]
Ni_range = [0.02, 0.12]
E_final_range = [0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]
R_range = [1000, 2000]
K_range = [0, 50, 100]
S_range = [1.0-sigma_S, 1.0+sigma_S]
T_range = [-10, 2]

ranges_dict = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range, 'S' : S_range, 'T' : T_range}


flat_sampler = pd.read_csv(flat_sampler_path, names=['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T'],
                           skiprows=(2*n_steps-1) * (n_walkers))


martinez_values_path = os.path.join('results', SN_name+'_martinez_values.csv')
martinez_values = list(pd.read_csv(martinez_values_path, names='val').T.iloc[0].astype(float))

martinez_fits_path = os.path.join('results', 'martinez_bestfit_models', SN_name+'.dat')
martinez_fits = pd.read_csv(martinez_fits_path, sep=r'\s+')
martinez_fits['VPH'] = martinez_fits['VPH'] / 10**5   # cm/s to km/s

data_lum_path = os.path.join('results', SN_name+'_martinez.csv')
# data_veloc_path = os.path.join('results', SN_name+'_expansion_velocities.csv')
data_lum = pd.read_csv(data_lum_path)
# data_veloc = pd.read_csv(data_veloc_path,
#                       usecols=['t_from_discovery', 'line', 'absorption_mean_velocity','absorption_std_velocity'])
# data_veloc.rename({'absorption_mean_velocity':'veloc', 'absorption_std_velocity':'dveloc'}, axis='columns', inplace=True)
# data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
# data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
# remove first point which seems like an artifact
# data_veloc = data_veloc.loc[data_veloc['t_from_discovery'] > 20]

# SN_data = {'lum': data_lum, 'veloc': data_veloc}
SN_data = {'lum': data_lum}


# last_rows = flat_sampler.iloc[-n_walkers:]

def get_param_results_dict(sampler_df):
    dict = {}
    for param in ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']:
        avg = np.average(sampler_df[param])
        sigma_lower, sigma_upper = np.percentile(sampler_df[param], [16, 84])
        dict[param] = avg
        dict[param + '_lower'] = avg - sigma_lower
        dict[param + '_upper'] = sigma_upper - avg
    return dict

def rounded_str(x):
    if not np.isinf(x) and np.abs(x) > 0.0000001:
        rounded = round(x, 2-int(np.floor(np.log10(abs(x)))))
        if rounded > 100:
            rounded = int(rounded)
    else:
        rounded = 0
    return str(rounded)


def result_text_from_dict(sampler_df):
    param_dict = get_param_results_dict(sampler_df)
    res_text = 'MCMC-SNEC fit\n\n'
    for param in ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']:
        if (param != 'K') & (param != 'R'):
            res_text += param + ': ' + rounded_str(param_dict[param]) + r'$\pm$ [' +\
                            rounded_str(param_dict[param+'_lower']) + ',' +\
                            rounded_str(param_dict[param+'_upper']) + ']\n'
    if (param_dict['K'] == 0) & (param_dict['R'] == 0):
        res_text += 'no CSM'
    else:
        res_text += 'K' + ': ' + rounded_str(param_dict['K']) + r'$\pm$ [' + \
                    rounded_str(param_dict['K_lower']) + ',' + \
                    rounded_str(param_dict['K_upper']) + ']\n'
        res_text += 'R' + ': ' + rounded_str(param_dict['R']) + r'$\pm$ [' + \
                    rounded_str(param_dict['R_lower']) + ',' + \
                    rounded_str(param_dict['R_upper']) + ']\n'
    return res_text

fit_plotting_times = np.concatenate((np.arange(1.0, 2.0, 0.1),
                                    np.arange(2.0, 10.0, 0.5),
                                    np.arange(10.0, 100.0, 2.0),
                                    np.arange(100.0, 150.0, 0.5),
                                    np.arange(150.0, 200.0, 2.0)))


def plot_lum_with_fit(data, sampler_df, ranges_dict, output_dir, n_walkers, SN_name, martinez_values):
    ranges_list = mcmc_snec.dict_to_list(ranges_dict)
    data_lum = data['lum']
    data_x = data_lum['t_from_discovery']
    data_x_moved = data_lum['t_from_discovery']
    data_y = data_lum['Lum']
    dy0 = data_lum['dLum0']
    dy1 = data_lum['dLum1']
    f_fit, ax = plt.subplots(figsize=(14, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    add_label = True
    # final fits from SNEC-MCMC
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        data_x_moved = fit_plotting_times - T
        y_fit = interp_lum.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
        if not isinstance(y_fit, str):
            # multiply whole graph by scaling factor
            y_fit = y_fit * S
            if add_label:
                ax.plot(fit_plotting_times, y_fit, alpha=0.2, color='purple', label='SNEC-MCMC fits')
                add_label = False
            else:
                ax.plot(fit_plotting_times, y_fit, alpha=0.3, color='purple')
    results_text = result_text_from_dict(sampler_df)
    ax.text(1.02, 0.83, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='purple')
    # martinez fits (according to SNEC)

    # martinez fits (according to SNEC)
    data_x_moved = fit_plotting_times - martinez_values[7]
    y_fit = interp_lum.snec_interpolator(martinez_values[0:6], ranges_list, data_x_moved)
    y_fit = y_fit * martinez_values[6]
    ax.plot(fit_plotting_times, y_fit, color='green', label='Martinez values on SNEC', linewidth=2.0)

    # martinez fits (according to their model)
    data_x_moved = martinez_fits['TIME'] + martinez_values[7]
    y_fit = np.power(10, martinez_fits['LOGL'])
    y_fit = y_fit * martinez_values[6]
    ax.plot(data_x_moved, y_fit, color='orange', label='Martinez fits', linewidth=2.0)

    results_text = 'Martinez results (by SNEC)'+'\n\n' + \
                   'Mzams: ' + str(martinez_values[0])+'\n' + \
                   'Ni: ' + str(martinez_values[1])+'\n' + \
                   'E: ' + str(martinez_values[2])+'\n' + \
                   'Mix: '+str(martinez_values[5])+'\n' + \
                   'S: ' + str(martinez_values[6]) + '\n' + \
                   'T: ' + str(martinez_values[7]) + '\n' + \
                   'noCSM'
    ax.text(1.02, 0.4, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='green')
    # real observations for the SN
    ax.errorbar(data_x, data_y, yerr=[dy0, dy1], marker='o', linestyle='None', color='k', label=SN_name+' observations')
    ax.set_xlim(-2, 200)
    ax.set_title(SN_name, fontsize=18)
    ax.legend()
    plt.tight_layout()
    # ax.set_yscale('linear')
    # f_fit.savefig(os.path.join(output_dir, str(result_t) + '_' + SN_name + '_lum.png'))
    # ax.set_ylim(float(0.5 * 10 ** 41), float(1.6 * 10 ** 43))
    ax.set_ylim(np.min(data_y) * 0.5, np.max(data_y) * 5)
    ax.set_yscale('log')
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, str(result_t) + '_' + SN_name + '_lum_log.png'))


def plot_veloc_with_fit(data, sampler_df, ranges_dict, output_dir, n_walkers, SN_name, martinez_values):
    ranges_list = mcmc_snec.dict_to_list(ranges_dict)
    data_veloc = data['veloc']
    data_x = data_veloc['t_from_discovery']
    data_y = data_veloc['veloc']
    dy = data_veloc['dveloc']
    f_fit, ax = plt.subplots(figsize=(14, 8))
    x_plotting = np.arange(0, 140, 0.5)
    add_label = True
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        data_x_moved = x_plotting - T
        y_fit = interp_veloc.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
        if not isinstance(y_fit, str):
            if add_label:
                ax.plot(x_plotting, y_fit, alpha=0.2, color='purple', label='SNEC-MCMC fits')
                add_label = False
            else:
                ax.plot(x_plotting, y_fit, alpha=0.3, color='purple')
    results_text = result_text_from_dict(sampler_df)
    ax.text(1.02, 0.83, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='purple')

    # martinez fits (according to SNEC)
    data_x_moved = x_plotting - martinez_values[7]
    y_fit = interp_veloc.snec_interpolator(martinez_values[0:6], ranges_list, data_x_moved)
    y_fit = y_fit * martinez_values[6]
    ax.plot(x_plotting, y_fit, color='green', label='Martinez values on SNEC', linewidth=2.0)

    # martinez fits (according to their model)
    data_x_moved = martinez_fits['TIME'] + martinez_values[7]
    ax.plot(data_x_moved, martinez_fits['VPH'], color='orange', label='Martinez fits', linewidth=2.0)

    results_text = 'Martinez results (by SNEC)' + '\n\n' + \
                   'Mzams: ' + str(martinez_values[0]) + '\n' + \
                   'Ni: ' + str(martinez_values[1]) + '\n' + \
                   'E: ' + str(martinez_values[2]) + '\n' + \
                   'Mix: ' + str(martinez_values[5]) + '\n' + \
                   'S: ' + str(martinez_values[6]) + '\n' + \
                   'T: ' + str(martinez_values[7]) + '\n' + \
                   'noCSM'
    ax.text(1.02, 0.4, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='green')
    # real observations for the SN
    ax.errorbar(data_x, data_y, yerr=dy, marker='o', linestyle='None', color='k',
                label=SN_name + ' observations')
    ax.legend()
    ax.set_xlim(-2, 140)
    ax.set_title(SN_name, fontsize=18)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, str(result_t) + '_' + SN_name + '_veloc.png'))


plot_lum_with_fit(SN_data, flat_sampler, ranges_dict,'figures', n_walkers, SN_name, martinez_values)

# plot_veloc_with_fit(SN_data, flat_sampler, ranges_dict,'figures', n_walkers, SN_name, martinez_values)





