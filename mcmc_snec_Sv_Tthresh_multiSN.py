import emcee
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import snec_result_interpolator_lum as interp_lum
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
T_thresh = 10 ** 3.75




def sigma_S(SN_name):
    distances = pd.read_csv(os.path.join('results', 'distances.csv'))
    distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
    distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
    sigma_S = 2 * distance_err / (distance)
    return sigma_S


def get_maximal_smaller(number, list):
    arr = np.array(list)
    smaller = arr[arr < number]
    return np.max(smaller)

def time_before_T_thresh(theta, ranges_dict):
    M_below = get_maximal_smaller(theta[0], ranges_dict['Mzams'])
    Ni_below = get_maximal_smaller(theta[1], ranges_dict['Ni'])
    E_below = get_maximal_smaller(theta[2], ranges_dict['E'])
    Mix_below = get_maximal_smaller(theta[5], ranges_dict['Mix'])
    model = 'M' + str(M_below) + '_Ni' + str(Ni_below) + '_E' + str(E_below) + \
            '_Mix' + str(Mix_below) + '_R0_K0'
    temps = pd.read_csv(os.path.join('..', 'all_temp_rad_data', model, 'T_eff.dat'), names=['time', 'temp'], sep=r'\s+')
    max_temp_below_Tthresh = np.max(temps['temp'].loc[temps['temp'] < T_thresh])
    time_thresh = float(temps['time'].loc[temps['temp'] == max_temp_below_Tthresh]) / 86400
    # time_thresh = 70
    return time_thresh

def theta_in_range(theta, ranges_dict):
    ranges_list = dict_to_list(ranges_dict)
    truth_value = True
    print(theta)
    print(ranges_dict)
    for i in range(len(ranges_list)):
        truth_value = truth_value and\
                      (ranges_list[i][0] <= theta[i] <= ranges_list[i][-1])
    return truth_value


def log_prior(theta, SN_names):
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
    num_SNe = len(SN_names)
    mu = 1.0
    log_prior = 0.0
    for i in range(num_SNe):
        S = theta[8*i+6]
        sigma = sigma_S(SN_names[i])
        log_prior += np.log(1.0 / (np.sqrt(2 * np.pi) * sigma)) - 0.5 * (S - mu) ** 2 / sigma ** 2
    # print('logprior', log_prior)
    return log_prior


def dict_to_list(dict):
    l = []
    for key in dict.keys():
        l.append(dict[key])
    return l


def calc_lum_likelihood(theta, data, ranges_dict, fitting_type):
    if theta_in_range(theta, ranges_dict):
        print('ok SN')
        ranges_list = dict_to_list(ranges_dict)
        if fitting_type == 'lum_veloc':
            num_obs = 1
        elif fitting_type == 'lum_veloc_normalized' \
            or fitting_type == 'lum_veloc_Tthresh_normalized':
            num_obs = len(data['t_from_discovery'])
        else:
            print('fitting_type should be: lum_veloc, lum_veloc_normalized or lum_veloc_normalized')
        data_x_moved = data['t_from_discovery'] - theta[7]
        data_y = data['Lum']
        data_dy = data['dLum0']
        y_fit = interp_lum.snec_interpolator(theta[0:6], ranges_list, data_x_moved)
        if not isinstance(y_fit, str):
            # multiply whole graph by scaling factor
            y_fit = y_fit * theta[6]
            # calculate the log likelihood
            log_likeli = - 0.5 * np.sum(np.log(2 * np.pi * data_dy ** 2) + ((data_y - y_fit) / data_dy) ** 2) / num_obs
        else:
            log_likeli = - np.inf
            print('impossible SN')
    else:
        log_likeli = - np.inf  # just a very big number so it won't go past the edge values
        print('out of range')
    return log_likeli

def calc_veloc_likelihood(theta, data, ranges_dict, fitting_type):
    if theta_in_range(theta, ranges_dict):
        print('ok SN')
        ranges_list = dict_to_list(ranges_dict)
        if fitting_type == 'lum_veloc':
            num_obs = 1
            data_veloc = data
        elif fitting_type == 'lum_veloc_normalized':
            num_obs = len(data['t_from_discovery'])
            data_veloc = data
        elif fitting_type == 'lum_veloc_Tthresh_normalized':
            max_veloc_time = time_before_T_thresh(theta, ranges_dict)
            data_veloc = data.loc[data['t_from_discovery'] <= max_veloc_time]
            num_obs = len(data_veloc['t_from_discovery'])
        else:
            print('fitting_type should be: lum_veloc, lum_veloc_normalized or lum_veloc_normalized')
            return np.nan
        if num_obs == 0:
            print('numobs is 0')
            return - np.inf
        else:
            data_x_moved = data_veloc['t_from_discovery'] - theta[7]
            data_y = data_veloc['veloc']
            data_dy = data_veloc['dveloc']
            y_fit = interp_veloc.snec_interpolator(theta[0:6], ranges_list, data_x_moved)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * theta[8]
                # calculate the log likelihood
                log_likeli = - 0.5 * np.sum(np.log(2 * np.pi * data_dy ** 2) + ((data_y - y_fit) / data_dy) ** 2) / num_obs
            else:
                log_likeli = - np.inf
                print('impossible SN')
    else:
        log_likeli = - np.inf  # just a very big number so it won't go past the edge values
        print('out of range')
    return log_likeli


def log_likelihood(theta, data, ranges_dict, fitting_type):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    ranges_dict : as provided by the user, in the form described above

    fitting_type : lum, veloc, combined or combined_normalized
    """
    log_likeli = 0
    num_SNe = len(ranges_dict)
    for i in range(num_SNe):
        SN_theta = np.concatenate((theta[8*i:8*(i+1)], [theta[num_SNe*8]]))
        print(SN_theta)
        log_likeli += calc_lum_likelihood(SN_theta, data[i]['lum'], ranges_dict[i], fitting_type) + \
                      calc_veloc_likelihood(SN_theta, data[i]['veloc'], ranges_dict[i], fitting_type)
    print('loglik', log_likeli)
    return log_likeli


def log_posterior(theta, data, ranges_dict, SN_names, fitting_type):
    lp = log_prior(theta, SN_names)
    ll = log_likelihood(theta, data, ranges_dict, fitting_type)
    log_post = lp + ll
    print('logpost', log_post)
    return log_post


def initial_guesses(ranges_dict, n_walkers):
    guesses = []
    num_SNe = len(ranges_dict)
    params = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
    for i in range(num_SNe):
        for param in params:
            guesses.append(np.random.rand(n_walkers) *
                       (ranges_dict[i][param][-1] - ranges_dict[i][param][0]) + ranges_dict[i][param][0])
    guesses.append(np.random.rand(n_walkers) *
                   (ranges_dict[0]['Sv'][-1] - ranges_dict[0]['Sv'][0]) + ranges_dict[0]['Sv'][0])
    return np.array(guesses).T


def emcee_fit_params(data, n_walkers, n_steps, ranges_dict, SN_names, fitting_type):
    n_params = len(SN_names) * 8 + 1
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=[data, ranges_dict, SN_names, fitting_type])
    init_guesses = initial_guesses(ranges_dict, n_walkers)
    print('init', init_guesses)
    sampler.run_mcmc(init_guesses, n_steps)
    return sampler


def corner_plot(sampler_chain_flat, ranges_dict, output_dir, SN_name):
    labels = list(ranges_dict.keys())
    corner_range = [1.] * len(labels)
    f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
    f_corner.savefig(os.path.join(output_dir, 'corner_plot'+SN_name+'.png'))


def chain_plots(sampler_chain, ranges_dict, output_dir, burn_in, SN_name):
    keys = list(ranges_dict.keys())
    for i in range(len(keys)):
        key = keys[i]
        plt.figure()
        plt.plot(sampler_chain[:, :, i].T, color='k', alpha=0.2)
        plt.xlabel('Step Number')
        plt.ylabel(key)
        plt.axvspan(0, burn_in, alpha=0.1, color='grey')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, key+'_'+SN_name+'.png'))


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


def get_param_results_dict(sampler_chain, ranges_dict, step, output_dir):
    params = list(ranges_dict.keys())
    dict = {}
    for i in range(len(params)):
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


def plot_lum_with_fit(data_lum, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, fitting_type):
    ranges_list = dict_to_list(ranges_dict)
    data_x = data_lum['t_from_discovery']
    data_y = data_lum['Lum']
    dy0 = data_lum['dLum0']
    dy1 = data_lum['dLum1']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    x_plotting = np.arange(0, 200, 0.5)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T, Sv] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T, Sv]
        if theta_in_range(requested, ranges_dict) or R == 0:
            data_x_moved = x_plotting - T
            y_fit = interp_lum.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                log_likeli_current = calc_lum_likelihood(requested, data_lum, ranges_dict, fitting_type)
                if log_likeli_current != - np.inf:
                    log_likeli.append(log_likeli_current)
                    ax.plot(x_plotting, np.log10(y_fit), alpha=0.4)
    log_likeli = np.average(log_likeli)
    data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 200)
    ax.set_ylim(40.8, 43)
    # ax.set_ylim(np.min(data_y) * 0.5, np.max(data_y) * 5)
    ax.set_title('step ' + str(step)
                 + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    # ax.set_ylim(float(3.0 * 10 ** 39), float(1.6 * 10 ** 43))
    # ax.set_yscale('log')
    # f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit' + str(step) + 'log.png'))
    # ax.set_yscale('linear')
    f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit_' + SN_name + '_' + str(step) + '.png'))
    return log_likeli


def plot_veloc_with_fit(data_veloc, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, fitting_type, Tthresh_normalized=False):
    ranges_list = dict_to_list(ranges_dict)
    data_x = data_veloc['t_from_discovery']
    data_y = data_veloc['veloc']
    dy = data_veloc['dveloc']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    log_likeli = []
    x_plotting = np.arange(0, np.max(data_veloc['t_from_discovery'])+5, 0.5)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T, Sv] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T, Sv]
        if theta_in_range(requested, ranges_dict) or R == 0:
            if Tthresh_normalized:
                max_veloc_time = time_before_T_thresh(requested, ranges_dict)
                thresh_veloc = data_veloc.loc[data_veloc['t_from_discovery'] <= max_veloc_time]
                x_plotting = np.arange(0, max_veloc_time+0.5, 0.5)
                data_x_moved = x_plotting - T
                y_fit = interp_veloc.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
                if not isinstance(y_fit, str):
                    # multiply whole graph by scaling factor
                    y_fit = y_fit * Sv
                    log_likeli_current = calc_veloc_likelihood(requested, thresh_veloc, ranges_dict, fitting_type)
                    if log_likeli_current != - np.inf:
                        log_likeli.append(log_likeli_current)
                        ax.plot(x_plotting, y_fit, alpha=0.4)
            else:
                data_x_moved = x_plotting - T
                y_fit = interp_veloc.snec_interpolator(requested[0:6], ranges_list, data_x_moved)
                if not isinstance(y_fit, str):
                    # multiply whole graph by scaling factor
                    y_fit = y_fit * Sv
                    log_likeli_current = calc_veloc_likelihood(requested, data_veloc, ranges_dict, fitting_type)
                    if log_likeli_current != - np.inf:
                        log_likeli.append(log_likeli_current)
                        ax.plot(x_plotting, y_fit, alpha=0.4)
    log_likeli = np.average(log_likeli)
    ax.errorbar(data_x, data_y, yerr=dy, marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 140)
    ax.set_title('step '+str(step)
                 +'\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, 'velocities_fit_' + SN_name + '_' + str(step) + '.png'))
    return log_likeli


def plot_lightcurve_with_fit(sampler_chain, SN_data, ranges_dict, fitting_type, output_dir, n_walkers, SN_name, step):
    if fitting_type == 'lum_veloc_Tthresh_normalized':
        log_likeli = 0
        log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, fitting_type)
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, fitting_type, Tthresh_normalized=True)
    else:
        log_likeli = 0
        log_likeli += plot_lum_with_fit(SN_data['lum'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, fitting_type)
        log_likeli += plot_veloc_with_fit(SN_data['veloc'], sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name, fitting_type)
    return log_likeli


