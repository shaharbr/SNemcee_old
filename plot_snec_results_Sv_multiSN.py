import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec_Sv_multiSN as mcmc_snec

import snec_result_interpolator_lum as interp_lum
import snec_result_interpolator_veloc as interp_veloc
from matplotlib import pyplot as plt




'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

SN_names = ['SN2004a', 'SN2005cs', 'SN2012aw', 'SN2012ec', 'SN2017eaw', 'SN2018aoq', 'SN2004et_d5.9']
low_luminosity_SNe = ['SN2005cs', 'SN2018aoq']


Mzams_range = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
Ni_range = [0.02, 0.07, 0.12]
E_final_range = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]
Mix_range = [2.0, 5.0, 8.0]

Mzams_range_low = [9.0, 10.0]
Ni_range_low = [0.001, 0.01]
E_final_range_low = [0.1, 0.3]
Mix_range_low = [5.0, 7.0]

R_range = [0, 0]
K_range = [0, 0]
S_range = [0.5, 2.0]
T_range = [-10, 2]
Sv_range = [0.1, 2.0]

param_high = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range, 'S' : S_range, 'T' : T_range,
                    'Sv' : Sv_range}

param_low = {'Mzams' : Mzams_range_low, 'Ni' : Ni_range_low, 'E' : E_final_range_low,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range_low, 'S' : S_range, 'T' : T_range,
                    'Sv' : Sv_range}

n_walkers = 116
n_steps = 30
n_params = 57
burn_in = 0

result_t = '2021-02-02_01-32-34_lum_veloc_findingSv_allSNe'
result_dir = os.path.join('mcmc_results', result_t)

flat_sampler_path = os.path.join(result_dir, 'flat_sampler.csv')
# flat_sampler = pd.read_csv(flat_sampler_path, header=None,
#                            skiprows=(n_steps-1) * (n_walkers))
flat_sampler_df = pd.read_csv(flat_sampler_path, header=None)
sampler = flat_sampler_df.to_numpy().reshape([n_walkers, n_steps, n_params])
flat_sampler = flat_sampler_df.to_numpy().reshape([n_steps * n_walkers, n_params])

data_lum_list = []
data_lum_w_early_list = []
data_veloc_list = []
parameter_ranges = []

for SN_name in SN_names:
    # import SN bolometric lum data
    data_lum_w_early = pd.read_csv(os.path.join('results', SN_name + '_martinez.csv'))
    # take only days between 30 and 200 (like the Martinez paper)
    data_lum_w_early = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] < 200]
    data_lum = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] > 30]
    # import SN photospheric velocities lum data
    data_veloc = pd.read_csv(os.path.join('results', SN_name + '_expansion_velocities.csv'),
                             usecols=['t_from_discovery', 'line', 'absorption_mean_velocity',
                                      'absorption_std_velocity'])
    data_veloc.rename({'absorption_mean_velocity': 'veloc', 'absorption_std_velocity': 'dveloc'}, axis='columns',
                      inplace=True)
    data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
    data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
    # remove first point which seems like an artifact
    data_veloc = data_veloc.loc[data_veloc['t_from_discovery'] > 20]
    # add to list of all SNe data
    data_lum_list.append(data_lum)
    data_lum_w_early_list.append(data_lum_w_early)
    data_veloc_list.append(data_veloc)
    sigma_S = mcmc_snec.sigma_S(SN_name)
    if SN_name in low_luminosity_SNe:
        param_low['S'] = [1.0 - sigma_S, 1.0 + sigma_S]
        parameter_ranges.append(param_low)
    else:
        param_high['S'] = [1.0 - sigma_S, 1.0 + sigma_S]
        parameter_ranges.append(param_high)


num_SNe = len(SN_names)

SN_data_all = [{'veloc': data_veloc_list[i], 'lum': data_lum_list[i]} for i in range(num_SNe)]
SN_data_all_w_early = [{'veloc': data_veloc_list[i], 'lum': data_lum_w_early_list[i]} for i in range(num_SNe)]

fit_plotting_times = np.concatenate((np.arange(1.0, 2.0, 0.1),
                                    np.arange(2.0, 10.0, 0.5),
                                    np.arange(10.0, 100.0, 2.0),
                                    np.arange(100.0, 150.0, 0.5),
                                    np.arange(150.0, 200.0, 2.0)))


# def get_param_results_dict(sampler_df, ranges_dict):
#     params = list(ranges_dict.keys())
#     dict = {}
#     for i in range(len(params)):
#         print(params[i])
#         last_results = sampler_df.iloc[:, i]
#         avg = np.average(last_results)
#         sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
#         dict[params[i]] = avg
#         dict[params[i]+'_lower'] = avg - sigma_lower
#         dict[params[i] + '_upper'] = sigma_upper - avg
#     return dict
#
#
# def rounded_str(x):
#     if not np.isinf(x) and np.abs(x) > 0.0000001:
#         rounded = round(x, 2-int(np.floor(np.log10(abs(x)))))
#         if rounded > 100:
#             rounded = int(rounded)
#     else:
#         rounded = 0
#     return str(rounded)
#
#
# def result_text_from_dict(sampler_df, ranges_dict):
#     param_dict = get_param_results_dict(sampler_df, ranges_dict)
#     res_text = 'MCMC-SNEC fit\n\n'
#     for param in ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']:
#         if (param != 'K') & (param != 'R'):
#             res_text += param + ': ' + rounded_str(param_dict[param]) + r'$\pm$ [' +\
#                             rounded_str(param_dict[param+'_lower']) + ',' +\
#                             rounded_str(param_dict[param+'_upper']) + ']\n'
#     if (param_dict['K'] == 0) & (param_dict['R'] == 0):
#         res_text += 'no CSM'
#     else:
#         res_text += 'K' + ': ' + rounded_str(param_dict['K']) + r'$\pm$ [' + \
#                     rounded_str(param_dict['K_lower']) + ',' + \
#                     rounded_str(param_dict['K_upper']) + ']\n'
#         res_text += 'R' + ': ' + rounded_str(param_dict['R']) + r'$\pm$ [' + \
#                     rounded_str(param_dict['R_lower']) + ',' + \
#                     rounded_str(param_dict['R_upper']) + ']\n'
#     return res_text
#
#
#
#
# def plot_lum_with_fit(data_lum, sampler_df, ranges_dict, output_dir, n_walkers, SN_name, fitting_type):
#     ranges_list = mcmc_snec.dict_to_list(ranges_dict)
#     data_x = data_lum['t_from_discovery']
#     data_y = data_lum['Lum']
#     dy0 = data_lum['dLum0']
#     dy1 = data_lum['dLum1']
#     f_fit, ax = plt.subplots(figsize=(10, 8))
#     ax.axvspan(-2, 30, alpha=0.1, color='grey')
#     log_likeli = []
#     x_plotting = np.arange(0, 200, 0.5)
#     print('start of walker loop')
#     for i in range(n_walkers):
#         [Mzams, Ni, E, R, K, Mix, S, T, Sv] = sampler_df.iloc[i]
#         requested = [Mzams, Ni, E, R, K, Mix, S, T, Sv]
#         print(requested)
#         print(ranges_list)
#         if mcmc_snec.theta_in_range(requested, ranges_dict) or R == 0:
#             data_x_moved = x_plotting - T
#             y_fit = interp_lum.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
#             if not isinstance(y_fit, str):
#                 multiply whole graph by scaling factor
                # y_fit = y_fit * S
                # ax.plot(x_plotting, np.log10(y_fit), alpha=0.4)
                # log_likeli.append(mcmc_snec.calc_lum_likelihood(requested, data_lum, ranges_dict, fitting_type))
    # log_likeli = np.average(log_likeli)
    # data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    # data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    # ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')
    # results_text = result_text_from_dict(sampler_df, ranges_dict)
    # ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
    #         verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    # ax.set_xlim(-2, 200)
    # ax.set_ylim(40.8, 43)
    # ax.set_title(SN_name + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    # plt.tight_layout()
    # f_fit.savefig(os.path.join(output_dir, str(result_t) + '_' + SN_name + '_lum_log.png'))
    # return log_likeli
#
#
# def plot_veloc_with_fit(data_veloc, sampler_df, ranges_dict, output_dir, n_walkers, SN_name, fitting_type, Tthresh_normalized=False):
#     ranges_list = mcmc_snec.dict_to_list(ranges_dict)
#     data_x = data_veloc['t_from_discovery']
#     data_y = data_veloc['veloc']
#     dy = data_veloc['dveloc']
#     f_fit, ax = plt.subplots(figsize=(10, 8))
#     log_likeli = []
#     x_plotting = np.arange(0, np.max(data_veloc['t_from_discovery'])+5, 0.5)
#     for i in range(n_walkers):
#         [Mzams, Ni, E, R, K, Mix, S, T, Sv] = sampler_df.iloc[i]
#         requested = [Mzams, Ni, E, R, K, Mix, S, T, Sv]
#         if mcmc_snec.theta_in_range(requested, ranges_dict) or R == 0:
#             if Tthresh_normalized:
#                 max_veloc_time = mcmc_snec.time_before_T_thresh(requested, ranges_dict)
#                 thresh_veloc = data_veloc.loc[data_veloc['t_from_discovery'] <= max_veloc_time]
#                 x_plotting = np.arange(0, max_veloc_time+0.5, 0.5)
#                 data_x_moved = x_plotting - T
#                 y_fit = interp_veloc.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
#                 if not isinstance(y_fit, str):
#                     multiply whole graph by scaling factor
                    # y_fit = y_fit * Sv
                    # ax.plot(x_plotting, y_fit, alpha=0.4)
                    # log_likeli.append(mcmc_snec.calc_veloc_likelihood(requested, thresh_veloc, ranges_dict, fitting_type))
            # else:
            #     data_x_moved = x_plotting - T
            #     y_fit = interp_veloc.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
            #     if not isinstance(y_fit, str):
            #         multiply whole graph by scaling factor
                    # y_fit = y_fit * Sv
                    # ax.plot(x_plotting, y_fit, alpha=0.4)
                    # log_likeli.append(mcmc_snec.calc_veloc_likelihood(requested, data_veloc, ranges_dict, fitting_type))
    # log_likeli = np.average(log_likeli)
    # ax.errorbar(data_x, data_y, yerr=dy, marker='o', linestyle='None', color='k')
    # results_text = result_text_from_dict(sampler_df, ranges_dict)
    # ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
    #         verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    # ax.set_xlim(-2, 140)
    # ax.set_title(SN_name + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    # plt.tight_layout()
    # f_fit.savefig(os.path.join(output_dir, str(result_t) + '_' + SN_name + '_veloc.png'))
    # return log_likeli
#
#
# def plot_lightcurve_with_fit(sampler_df, SN_data, ranges_dict, fitting_type, output_dir, n_walkers, SN_name):
#     if fitting_type == 'lum_veloc_Tthresh_normalized':
#         log_likeli = 0
#         log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_df, ranges_dict, output_dir, n_walkers, SN_name, fitting_type)
#         log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_df, ranges_dict, output_dir, n_walkers, SN_name, fitting_type, Tthresh_normalized=True)
#     else:
#         log_likeli = 0
#         log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_df, ranges_dict, output_dir, n_walkers, SN_name, fitting_type)
#         log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_df, ranges_dict, output_dir, n_walkers, SN_name, fitting_type)
#     return log_likeli
#
# def corner_plot(sampler_df, ranges_dict, output_dir, SN_name):
#     labels = list(ranges_dict.keys())
#     corner_range = [1.] * len(labels)
#     f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
#     f_corner.savefig(os.path.join(output_dir, 'corner_plot'+SN_name+'.png'))
#
# def chain_plots(sampler_df, ranges_dict, output_dir, burn_in, SN_name):
#     keys = list(ranges_dict.keys())
#     for i in range(len(keys)):
#         key = keys[i]
#         plt.figure()
#         plt.plot(sampler_chain[:, :, i].T)
#         plt.xlabel('Step Number')
#         plt.ylabel(key)
#         plt.axvspan(0, burn_in, alpha=0.1, color='grey')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir, key+'_'+SN_name+'.png'))

res_dir = 'figures'


for i in range(num_SNe):
    sampler_8p = sampler[:, :, 8 * i:8 * (i + 1)]
    sampler_1p = sampler[:, :, 56].reshape(n_walkers,n_steps, 1)
    concat_sampler = np.concatenate((sampler_8p, sampler_1p), axis=2)

    flat_sampler_8p = flat_sampler[n_walkers * burn_in:-1, 8 * i:8 * (i + 1)]
    flat_sampler_1p = flat_sampler[n_walkers * burn_in:-1, 56].reshape(n_walkers * (n_steps - burn_in) - 1, 1)
    concat_flat_sampler = np.concatenate((flat_sampler_8p, flat_sampler_1p), axis=1)

    mcmc_snec.plot_lightcurve_with_fit(concat_sampler,
        SN_data_all_w_early[i], parameter_ranges[i],
        'lum_veloc_normalized', res_dir, n_walkers, SN_names[i], n_steps - 1)

    mcmc_snec.plot_lightcurve_with_fit(concat_sampler,
        SN_data_all_w_early[i], parameter_ranges[i],
        'lum_veloc_normalized', res_dir, n_walkers, SN_names[i], 0)
    mcmc_snec.chain_plots(concat_sampler,
                          parameter_ranges[i], res_dir, burn_in, SN_names[i])

    mcmc_snec.corner_plot(concat_flat_sampler,
                          parameter_ranges[i], res_dir, SN_names[i])
