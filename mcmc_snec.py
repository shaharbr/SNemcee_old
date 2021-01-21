import emcee
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import snec_result_interpolator_lum as interp_lum
import snec_result_interpolator_mag as interp_mag
import snec_result_interpolator_veloc as interp_veloc
import pandas as pd
import os
import corner
import csv
import copy


'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'

colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}

# multiplicative factor for SNEC's predicted photospheric velocities to Fe II velocities,
# found to be equal about 1.4, by fitting to all the SNe in our sample
Sv = 1.4

def theta_in_range(theta, ranges_dict):
    ranges_list = dict_to_list(ranges_dict)
    truth_value = True
    for i in range(len(ranges_list)):
        truth_value = truth_value and\
                      (ranges_list[i][0] <= theta[i] <= ranges_list[i][-1])
    return truth_value


def log_prior(theta, ranges_dict, nonuniform_priors=None):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    ranges_dict : as provided by the user, in the form of (example):
    ranges_dict =  {Mzams : [13.0, 14.0, 15.0, 16.0],
                    Ni : [0.02, 0.07, 0.12, 0.17],
                    E : [0.7, 0.9, 1.1, 1.3, 1.5],
                    Mix : [5.0, 8.0],
                    R : [0.1, 0.1],
                    K : [0.1, 0.1],
                    S : [0.8, 1.2],
                    T : [-15, +15]}
    """
    if not theta_in_range(theta, ranges_dict):
        return -np.inf
    else:
        if nonuniform_priors is None:
            # flat probability for all means p=1 for all, and log_p=0, so sum is also 0.
            return 0.0
        else:
            param_dict = {'Mzams': theta[0], 'Ni': theta[1], 'E': theta[2], 'R':theta[3],
                          'K': theta[4], 'Mix':theta[5], 'S':theta[6], 'T':theta[7]}
            log_prior = 0.0
            for param in list(nonuniform_priors.keys()):
                if nonuniform_priors[param] == 'gaussian':
                    mu = nonuniform_priors[param]['gaussian']['mu']
                    sigma = nonuniform_priors[param]['gaussian']['sigma']
                    log_prior += np.log(1.0 / (np.sqrt(2 * np.pi) * sigma)) - 0.5 * (param_dict[param] - mu) ** 2 / sigma ** 2
                    print('logprior', log_prior)
                if nonuniform_priors[param] == 'polynomial':
                    log_prior += np.log(nonuniform_priors[param]['polynomial'](param_dict[param]))
            return log_prior


def dict_to_list(dict):
    l = []
    for key in dict.keys():
        l.append(dict[key])
    return l


def calc_lum_likelihood(theta, data, ranges_list):
    data_x_moved = data['t_from_discovery'] - theta[7]
    data_y = data['Lum']
    data_dy = data['dLum0']
    y_fit = interp_lum.snec_interpolator(theta[0:6], ranges_list, data_x_moved)
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]
        # calculate the log likelihood
        log_likeli = - 0.5 * np.sum(np.log(2 * np.pi * data_dy ** 2) + ((data_y - y_fit) / data_dy) ** 2)
    else:
        log_likeli = - np.inf
        print('impossible SN')
    return log_likeli

def calc_veloc_likelihood(theta, data, ranges_list):

    data_x_moved = data['t_from_discovery'] - theta[7]
    data_y = data['veloc']
    data_dy = data['dveloc']
    y_fit = interp_veloc.snec_interpolator(theta[0:6], ranges_list, data_x_moved)
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * Sv
        # calculate the log likelihood
        log_likeli = - 0.5 * np.sum(np.log(2 * np.pi * data_dy ** 2) + ((data_y - y_fit) / data_dy) ** 2)
    else:
        log_likeli = - np.inf
        print('impossible SN')
    return log_likeli


def calc_mag_likelihood(theta, data, ranges_list):
    log_likeli = 0
    data_x_allfilt = data['t_from_discovery'] - theta[7]
    filters = list(data['filter'].unique())
    y_fit = interp_mag.snec_interpolator(theta[0:6], ranges_list, data_x_allfilt, filters)
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]
        y_fit['time'] = data_x_allfilt
        for filt in filters:
            data_filt = data.loc[data['filter'] == filt]
            data_x_filt = data_filt['t_from_discovery'] - theta[7]
            data_y_filt = data_filt['abs_mag']
            data_dy_filt = data_filt['dmag']
            y_fit_filt = y_fit[filt].loc[y_fit['time'].isin(data_x_filt)]
            if not len(data_x_filt) == len(y_fit_filt):
                print('stuck')
            # calculate the log likelihood
            log_likeli += - 0.5 * np.sum(
                np.log(2 * np.pi * data_dy_filt ** 2) + ((data_y_filt - y_fit_filt) / data_dy_filt) ** 2)
            # print('log likelihood', log_likeli)
    else:
        log_likeli = - np.inf
        print('impossible SN')
    return log_likeli



def log_likelihood(theta, data, ranges_dict, fitting_type):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    ranges_dict : as provided by the user, in the form described above

    fitting_type : lum, mag, veloc, combined or combined_normalized
    """
    print(theta)
    ranges_list = dict_to_list(ranges_dict)
    print(ranges_list)
    if theta_in_range(theta, ranges_dict) or theta[3] == 0:
        print('ok SN')
        if fitting_type == 'lum':
            log_likeli = calc_lum_likelihood(theta, data['lum'], ranges_list)
        elif fitting_type == 'mag':
            log_likeli = calc_mag_likelihood(theta, data['mag'], ranges_list)
        elif fitting_type == 'veloc':
            log_likeli = calc_veloc_likelihood(theta, data['veloc'], ranges_list)
        elif fitting_type == 'lum_veloc':
            log_likeli = calc_lum_likelihood(theta, data['lum'], ranges_list) + \
                         calc_veloc_likelihood(theta, data['veloc'], ranges_list)
        elif fitting_type == 'lum_veloc_normalized':
            num_obs_lum = len(data['lum']['t_from_discovery'])
            num_obs_veloc = len(data['veloc']['t_from_discovery'])
            log_likeli = calc_lum_likelihood(theta, data['lum'], ranges_list) / num_obs_lum+ \
                         calc_veloc_likelihood(theta, data['veloc'], ranges_list) / num_obs_veloc
        elif fitting_type == 'mag_veloc':
            log_likeli = calc_mag_likelihood(theta, data['mag'], ranges_list) + \
                         calc_veloc_likelihood(theta, data['veloc'], ranges_list)
        elif fitting_type == 'mag_veloc_normalized':
            num_obs_mag = len(data['mag']['t_from_discovery'])
            num_obs_veloc = len(data['veloc']['t_from_discovery'])
            log_likeli = calc_mag_likelihood(theta, data['mag'], ranges_list) / num_obs_mag+ \
                         calc_veloc_likelihood(theta, data['veloc'], ranges_list) / num_obs_veloc
        elif fitting_type == 'combined':
            log_likeli = calc_lum_likelihood(theta, data['lum'], ranges_list) + \
                         calc_mag_likelihood(theta, data['mag'], ranges_list) + \
                         calc_veloc_likelihood(theta, data['veloc'], ranges_list)
        elif fitting_type == 'combined_normalized':
            num_obs_lum = len(data['lum']['t_from_discovery'])
            num_obs_mag = len(data['mag']['t_from_discovery'])
            num_obs_veloc = len(data['veloc']['t_from_discovery'])
            log_likeli = calc_lum_likelihood(theta, data['lum'], ranges_list) / num_obs_lum+ \
                         calc_mag_likelihood(theta, data['mag'], ranges_list) / num_obs_mag+ \
                         calc_veloc_likelihood(theta, data['veloc'], ranges_list) / num_obs_veloc
        else:
            print('fitting_type should be: lum, mag, veloc, lum_veloc, lum_veloc_normalized, mag_veloc, mag_veloc_normalized, combined or combined_normalized')
    else:
        log_likeli = - np.inf  # just a very big number so it won't go past the edge values
        print('out of range')
    print('loglik', log_likeli)
    return log_likeli



def log_posterior(theta, data, ranges_dict, fitting_type, nonuniform_priors=None):
    if theta_in_range(theta, ranges_dict):
        lp = log_prior(theta, ranges_dict, nonuniform_priors)
        ll = log_likelihood(theta, data, ranges_dict, fitting_type)
        log_post = lp + ll
    else:
        log_post = -np.inf
    print('logpost', log_post)
    return log_post


def initial_guesses(ranges_dict, n_walkers):
    guesses = []
    for param in list(ranges_dict.keys()):
        guesses.append(np.random.rand(n_walkers) *
                       (ranges_dict[param][-1] - ranges_dict[param][0]) + ranges_dict[param][0])
    return np.array(guesses).T


def emcee_fit_params(data, n_walkers, n_steps, ranges_dict, fitting_type, nonuniform_priors=None, init_guesses=None):
    n_params = len(ranges_dict.keys())
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=[data, ranges_dict, fitting_type, nonuniform_priors])
    if init_guesses is None:
        init_guesses = initial_guesses(ranges_dict, n_walkers)
    sampler.run_mcmc(init_guesses, n_steps)
    return sampler


def corner_plot(sampler_chain_flat, ranges_dict, output_dir):
    labels = list(ranges_dict.keys())
    corner_range = [1., 1., 1., 1., 1., 1., 1., 1.]
    print('shapeshapecorner', sampler_chain_flat.shape)
    f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
    f_corner.savefig(os.path.join(output_dir, 'corner_plot.png'))


def chain_plots(sampler_chain, ranges_dict, output_dir, burn_in, first_stage_steps=None):
    keys = list(ranges_dict.keys())
    for i in range(len(keys)):
        key = keys[i]
        plt.figure()
        plt.plot(sampler_chain[:, :, i].T)
        plt.xlabel('Step Number')
        plt.ylabel(key)
        plt.axvspan(0, burn_in, alpha=0.1, color='grey')
        if first_stage_steps is not None:
            plt.axvline(x=first_stage_steps, color='black')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, key+'.png'))


def get_each_walker_result(sampler_chain, ranges_dict, step):
    params = list(ranges_dict.keys())
    dict = {}
    for i in range(len(params)):
        last_results = sampler_chain[:, step:, i]
        avg = np.average(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg
    print(dict)

# def param_results_to_gaussians(results_dict):
#     dict = {}
#     params = list(results_dict.keys())
#     for param in params:
#         dict3





def get_param_results_dict(sampler_chain, ranges_dict, step, output_dir):
    params = list(ranges_dict.keys())
    dict = {}
    for i in range(len(params)):
        print(params[i])
        last_results = sampler_chain[:, step:, i]
        avg = np.average(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        # dict[params[i]+'_all_walkers'] = last_results
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg
    print(dict)
    with open(os.path.join(output_dir, 'final_results.csv'), 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)
    return dict

def rounded_str(x):
    if not np.isinf(x) and np.abs(x) > 0.0000001:
        rounded = round(x, 2-int(np.floor(np.log10(abs(x)))))
        if rounded > 100:
            rounded = int(rounded)
    else:
        rounded = 0
    return str(rounded)


def result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir):
    param_dict = get_param_results_dict(sampler_chain, ranges_dict, step, output_dir)
    res_text = SN_name + '\n'
    params = list(ranges_dict.keys())
    for param in params:
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


def plot_lum_with_fit(data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name):
    ranges_list = dict_to_list(ranges_dict)
    data_lum = data['lum']
    data_x = data_lum['t_from_discovery']
    data_y = data_lum['Lum']
    dy0 = data_lum['dLum0']
    dy1 = data_lum['dLum1']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    x_plotting = np.arange(0, 200, 0.5)
    print('start of walker loop')
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        print(requested)
        print(ranges_list)
        if theta_in_range(requested, ranges_dict) or R == 0:
            data_x_moved = x_plotting - T
            y_fit = interp_lum.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                ax.plot(x_plotting, y_fit, alpha=0.4)
                log_likeli.append(log_likelihood(requested, data, ranges_dict, 'lum'))
                print(log_likelihood(requested, data, ranges_dict, 'lum'))
    log_likeli = np.average(log_likeli)
    ax.errorbar(data_x, data_y, yerr=[dy0, dy1], marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 200)
    ax.set_ylim(np.min(data_y) * 0.5, np.max(data_y) * 5)
    ax.set_title('step ' + str(step)
                 + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    # ax.set_ylim(float(3.0 * 10 ** 39), float(1.6 * 10 ** 43))
    ax.set_yscale('log')
    f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit' + str(step) + 'log.png'))
    ax.set_yscale('linear')
    f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit' + str(step) + '.png'))
    return log_likeli

def plot_mag_with_fit(data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name):
    ranges_list = dict_to_list(ranges_dict)
    data_mag = data['mag']
    data_x = data_mag['t_from_discovery']
    filters = list(data_mag['filter'].unique())
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    x_plotting = np.arange(0, 200, 0.5)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict) or R == 0:
            data_x_moved_all = x_plotting - T
            y_fit = interp_mag.snec_interpolator(requested[0:6], ranges_list, data_x_moved_all, filters)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                for filt in filters:
                    ax.plot(x_plotting, y_fit[filt], color=colors[filt], alpha=0.3)
                log_likeli.append(log_likelihood(requested, data, ranges_dict, 'mag'))
    log_likeli = np.average(log_likeli)
    for filt in filters:
        data_filt = data_mag.loc[data_mag['filter'] == filt]
        data_x = data_filt['t_from_discovery']
        data_y = data_filt['abs_mag']
        data_dy = data_filt['dmag']
        ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None',
                    label=filt, color=colors[filt])
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_title('step ' + str(step)
                 + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    ax.set_xlim(-2, 137)
    ax.set_ylim(14, 20)
    f_fit.savefig(os.path.join(output_dir, 'mag_lightcurve_fit' + str(step) + '.png'))
    return log_likeli



def plot_veloc_with_fit(data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name):
    ranges_list = dict_to_list(ranges_dict)
    data_veloc = data['veloc']
    data_x = data_veloc['t_from_discovery']
    data_y = data_veloc['veloc']
    dy = data_veloc['dveloc']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    log_likeli = []
    x_plotting = np.arange(0, 200, 0.5)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict) or R == 0:
            data_x_moved = x_plotting - T
            y_fit = interp_veloc.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * Sv
                ax.plot(x_plotting, y_fit, alpha=0.4)
                log_likeli.append(log_likelihood(requested, data, ranges_dict, 'veloc'))
    log_likeli = np.average(log_likeli)
    ax.errorbar(data_x, data_y, yerr=dy, marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 140)
    ax.set_title('step '+str(step)
                 +'\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, 'velocities_fit' +str(step) +'.png'))
    return log_likeli



def plot_lightcurve_with_fit(sampler_chain, SN_data, ranges_dict, fitting_type, output_dir, n_walkers, SN_name, step):
    if fitting_type == 'lum':
        log_likeli = plot_lum_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'mag':
        log_likeli = plot_mag_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'veloc':
        log_likeli = plot_veloc_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'lum_veloc':
        log_likeli = 0
        log_likeli += plot_lum_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
        log_likeli += plot_veloc_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'lum_veloc_normalized':
        log_likeli = 0
        num_obs_veloc = len(SN_data['veloc']['t_from_discovery'])
        num_obs_lum = len(SN_data['lum']['t_from_discovery'])
        log_likeli += plot_lum_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_lum
        log_likeli += plot_veloc_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_veloc
    if fitting_type == 'mag_veloc':
        log_likeli = 0
        log_likeli += plot_mag_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
        log_likeli += plot_veloc_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'mag_veloc_normalized':
        log_likeli = 0
        num_obs_veloc = len(SN_data['veloc']['t_from_discovery'])
        num_obs_mag = len(SN_data['mag']['t_from_discovery'])
        log_likeli += plot_mag_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_mag
        log_likeli += plot_veloc_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_veloc
    if fitting_type == 'combined':
        log_likeli = 0
        log_likeli += plot_lum_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
        log_likeli += plot_mag_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
        log_likeli += plot_veloc_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name)
    if fitting_type == 'combined_normalized':
        log_likeli = 0
        num_obs_mag = len(SN_data['mag']['t_from_discovery'])
        num_obs_veloc = len(SN_data['veloc']['t_from_discovery'])
        num_obs_lum = len(SN_data['lum']['t_from_discovery'])
        log_likeli += plot_mag_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_mag
        log_likeli += plot_lum_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_lum
        log_likeli += plot_veloc_with_fit(SN_data, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name) / num_obs_veloc
    return log_likeli



