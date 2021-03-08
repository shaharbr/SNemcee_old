import emcee
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import snec_model_interpolator as interp
import pandas as pd
import os
import corner
import csv
from scipy.stats import chi2


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

# initializing global variable 'models' (a dictionary maintaining in cache memory all the loaded models)
models = {}


def emcee_fit_params(data_lum, n_walkers, n_steps, ranges_dict, nonuniform_priors=None, init_guesses=None):
    data_no_early_lum = data_lum.loc[data_lum['t_from_discovery'] > 30]
    n_params = len(ranges_dict.keys())
    if init_guesses is None:
        init_guesses = initial_guesses(ranges_dict, n_walkers)
    initialize_empty_models_dict(ranges_dict)
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior,
                                    args=[data_no_early_lum, ranges_dict, nonuniform_priors])
    sampler.run_mcmc(init_guesses, n_steps)
    return sampler


def log_prior(theta, ranges_dict, nonuniform_priors=None):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T, Sv].

    ranges_dict : dictionary with the theta parameters as keys, and lists of the parameter's sampled
                    values as as dictionary values.
                    Example:
                    example_dict =  {Mzams : [13.0, 14.0, 15.0, 16.0],
                                    Ni : [0.02, 0.07, 0.12, 0.17],
                                    E : [0.7, 0.9, 1.1, 1.3, 1.5],
                                    Mix : [5.0, 8.0],
                                    R : [0.1, 0.1],
                                    K : [0.1, 0.1],
                                    S : [0.8, 1.2],
                                    T : [-15, +15]}

    nonuniform_priors: a dictionary with the parameters for the gaussian or polynomial curve
                        describing the prior for each of the parameters in the dictionary keys.
                        Example:
                        Example_nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}},
                                                    'Mzams': {'polynomial': poly1d_object}}
    """
    if not theta_in_range(theta, ranges_dict):
        return -np.inf
    else:
        if nonuniform_priors is None:
            # flat probability for all means p=1 for all, and log_p=0, so sum is also 0.
            return 0.0
        else:
            param_dict = {'Mzams': theta[0], 'Ni': theta[1], 'E': theta[2], 'R': theta[3],
                          'K': theta[4], 'Mix': theta[5], 'S': theta[6], 'T': theta[7]}
            log_prior = 0.0
            for param in list(nonuniform_priors.keys()):
                if nonuniform_priors[param] == 'gaussian':
                    mu = nonuniform_priors[param]['gaussian']['mu']
                    sigma = nonuniform_priors[param]['gaussian']['sigma']
                    log_prior += np.log(1.0 / (np.sqrt(2 * np.pi) * sigma)) - 0.5 * (
                                param_dict[param] - mu) ** 2 / sigma ** 2
                    print('logprior', log_prior)
                if nonuniform_priors[param] == 'polynomial':
                    log_prior += np.log(nonuniform_priors[param]['polynomial'](param_dict[param]))
            return log_prior


def log_likelihood(theta, data_lum, ranges_dict):
    """
    Parameters
    ----------
    theta: vector containing a specific set of parameters theta = [Mzams, Ni, E, R, K, Mix, S, T].

    data_lum: luminosity observations DataFrame, with the headers: "t_from_discovery", "Lum", "dLum0", "dLum1"

    ranges_dict : as provided by the user, in the form described above

    """
    print(theta)
    if theta_in_range(theta, ranges_dict) or theta[3] == 0:
        print('ok SN')
        surrounding_values = get_surrouding_values(theta[0:6], ranges_dict)
        load_surrounding_models(theta[0:6], ranges_dict)
        log_likeli = calc_lum_likelihood(theta, data_lum, surrounding_values)
        print('loglik', log_likeli)
        return log_likeli
    else:
        print('out of range')
        return - np.inf  # just a very big number so it won't go past the edge values



def calc_lum_likelihood(theta, data_lum, surrounding_values):
    data_x_moved = data_lum['t_from_discovery'] - theta[7]
    data_y = data_lum['Lum']
    data_dy = data_lum['dLum0']
    y_fit = interp.snec_interpolator(theta[0:6], surrounding_values, models['lum'], data_x_moved)
    if not isinstance(y_fit, str):
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]
        # calculate the log likelihood
        df = len(data_y) - 1
        log_likeli = chi2.logpdf(chi_square_weighted(data_y, data_dy, y_fit), df)
        print(log_likeli)
        return log_likeli
    else:
        print('impossible SN')
        return - np.inf


def log_posterior(theta, data_lum, ranges_dict, nonuniform_priors=None):
    if theta_in_range(theta, ranges_dict):
        lp = log_prior(theta, ranges_dict, nonuniform_priors)
        ll = log_likelihood(theta, data_lum, ranges_dict)
        log_post = lp + ll
        print('logpost', log_post)
        return log_post
    else:
        print('out of range')
        return - np.inf  # just a very big number so it won't go past the edge values


def initialize_empty_models_dict(ranges_dict):
    models['lum'] = {}
    for Mzams in ranges_dict['Mzams']:
        models['lum'][Mzams] = {}
        for Ni in ranges_dict['Ni']:
            models['lum'][Mzams][Ni] = {}
            for E in ranges_dict['E']:
                models['lum'][Mzams][Ni][E] = {}
                for R in ranges_dict['R']:
                    models['lum'][Mzams][Ni][E][R] = {}
                    for K in ranges_dict['K']:
                        models['lum'][Mzams][Ni][E][R][K] = {}
                        for Mix in ranges_dict['Mix']:
                            models['lum'][Mzams][Ni][E][R][K][Mix] = None


def get_surrouding_values(requested, ranges_dict):
    params = list(ranges_dict.keys())
    surrouding_values = {param: [] for param in params}
    for i in range(len(requested)):
        param_range = np.array(ranges_dict[params[i]])
        below = np.max(param_range[param_range <= requested[i]])
        above = np.min(param_range[param_range >= requested[i]])
        surrouding_values[params[i]] = [below, above]
    return surrouding_values


def load_model(Mzams, Ni, E, R, K, Mix):
    if R == 0 or K == 0:
        R = 0
        K = 0
    name = 'M' + str(Mzams) + \
           '_Ni' + str(Ni) + \
           '_E' + str(E) + \
           '_Mix' + str(Mix) + \
           '_R' + str(R) + \
           '_K' + str(K)
    modelpath = os.path.join('..', 'all_lum_data', name, 'lum_observed.dat')
    if os.stat(modelpath).st_size < 10 ** 5:
        return 'failed SN'
    else:
        snec_model = pd.read_csv(modelpath,
                                 names=['t_from_discovery', 'Lum'], sep=r'\s+')
        time_col = snec_model['t_from_discovery'] / 86400  # sec to days
        interp_days = np.linspace(0, 196, 1961)
        snec_model = np.interp(interp_days, time_col, snec_model['Lum'])
        return snec_model


def load_surrounding_models(requested, ranges_dict):
    surrouding_values = get_surrouding_values(requested, ranges_dict)
    for Mzams in surrouding_values['Mzams']:
        for Ni in surrouding_values['Ni']:
            for E in surrouding_values['E']:
                for R in surrouding_values['R']:
                    for K in surrouding_values['K']:
                        for Mix in surrouding_values['Mix']:
                            if models['lum'][Mzams][Ni][E][R][K][Mix] is None:
                                models['lum'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix)


def theta_in_range(theta, ranges_dict):
    ranges_list = dict_to_list(ranges_dict)
    truth_value = True
    for i in range(len(ranges_list)):
        truth_value = truth_value and\
                      (ranges_list[i][0] <= theta[i] <= ranges_list[i][-1])
    return truth_value


def dict_to_list(dict):
    l = []
    for key in dict.keys():
        l.append(dict[key])
    return l


def chi_square_weighted(y, dy, y_fit):
    return np.sum(((y - y_fit) / dy) ** 2)


def initial_guesses(ranges_dict, n_walkers):
    guesses = []
    for param in list(ranges_dict.keys()):
        guesses.append(np.random.rand(n_walkers) *
                       (ranges_dict[param][-1] - ranges_dict[param][0]) + ranges_dict[param][0])
    return np.array(guesses).T


def corner_plot(sampler_chain_flat, ranges_dict, output_dir):
    labels = list(ranges_dict.keys())
    corner_range = [1.] * len(labels)
    f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
    f_corner.savefig(os.path.join(output_dir, 'corner_plot.png'))


def chain_plots(sampler_chain, ranges_dict, output_dir, burn_in, first_stage_steps=None):
    keys = list(ranges_dict.keys())
    for i in range(len(keys)):
        key = keys[i]
        plt.figure()
        plt.plot(sampler_chain[:, :, i].T, color='k', alpha=0.2)
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


def plot_lum_with_fit(data_lum, sampler_chain, ranges_dict, step, output_dir, n_walkers, SN_name):
    data_x = data_lum['t_from_discovery']
    data_y = data_lum['Lum']
    dy0 = data_lum['dLum0']
    dy1 = data_lum['dLum1']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    x_plotting = np.linspace(0, 196, 1961)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if theta_in_range(requested, ranges_dict) or R == 0:
            data_x_moved = x_plotting - T
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['lum'], data_x_moved)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                ax.plot(x_plotting, np.log10(y_fit), alpha=0.4)
                log_likeli.append(calc_lum_likelihood(requested, data_lum, surrounding_values))
    log_likeli = np.average(log_likeli)
    data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')
    results_text = result_text_from_dict(sampler_chain, ranges_dict, SN_name, step, output_dir)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(-2, 196)
    ax.set_ylim(40.8, 43)
    ax.set_title('step ' + str(step)
                 + '\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, 'lum_lightcurve_fit' + str(step) + '.png'))
    return log_likeli



