import emcee
import numpy as np
from matplotlib import pyplot as plt
import snec_result_interpolator_mag as interp
import pandas as pd
import datetime
from pathlib import Path
import os
import corner
import csv


'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

SN_name = 'SN2004et'
Mzams_range = [13.0, 14.0, 15.0, 16.0]
Ni_range = [0.02, 0.07, 0.12, 0.17]
E_final_range = [0.7, 0.9, 1.1, 1.3, 1.5]
Mix_range = [5.0, 8.0]
R_range = [0.1, 0.1]
K_range = [0.1, 0.1]
S_range = [0.8, 1.2]
T_range = [0, 30] # because can't have negative values, do 15 minus diff (so 0 is -15, and 30 is +15)
pysynphot_models = True

n_walkers = 20
n_steps = 10
n_params = 8
burn_in = 0


time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_mag')
Path(res_dir).mkdir(parents=True, exist_ok=True)

run_param_df = pd.DataFrame.from_dict({'Mzams_range': str(Mzams_range), 'Ni_range': str(Ni_range),
                             'E_final_range': str(E_final_range), 'R_range': str(R_range), 'K_range': str(K_range),
                             'Mix_range': str(Mix_range), 'Scaling_range': str(S_range), 'T_exp_min15_range': str(T_range),
                             'n_walkers': n_walkers, 'n_steps': n_steps, 'n_params': n_params,
                             'burn_in': burn_in,
                             'time': time_now,
                             'pysynphot_model': str(pysynphot_models)}, orient='index')
run_param_df.to_csv(os.path.join(res_dir, 'run_parameters.csv'))

m_Solar = 1.989 * (10 ** 33)  # gram

# import SN mag data
data_filepath = os.path.join('results', SN_name+'_lightcurves')
SN_w_early = pd.read_csv(data_filepath, usecols=['dmag', 'filter', 'abs_mag', 't_from_discovery'])
if pysynphot_models:
    SN_w_early = SN_w_early.loc[SN_w_early['t_from_discovery'] < 120]
SN_w_early = SN_w_early.loc[SN_w_early['t_from_discovery'] < 200]

SN_w_early['abs_mag'] = SN_w_early['abs_mag'].abs()
filters = list(SN_w_early['filter'].unique())
filters = list(set(filters).intersection(['g', 'r', 'i', 'V', 'R', 'I']))
SN_w_early = SN_w_early.loc[SN_w_early['filter'].isin(filters)]
SN_w_early = SN_w_early.sort_values('t_from_discovery')
SN = SN_w_early.loc[SN_w_early['t_from_discovery'] > 30]

colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}


plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'


def log_prior(theta):
    # theta is a vector containing a specific set of parameters theta = [M, Ni, E, R, K, Mix, S, T]

    # Mzams (M)
    if Mzams_range[0] <= theta[0] <= Mzams_range[-1]:
        prob_M = 1. / theta[0]
    else:
        prob_M = 0.

    # Ni mass (Ni)
    if Ni_range[0] <= theta[1] <= Ni_range[-1]:
        prob_Ni = 1. / theta[1]
    else:
        prob_Ni = 0.

    # E_final (E)
    if E_final_range[0] <= theta[2] <= E_final_range[-1]:
        prob_E = 1. / theta[2]
    else:
        prob_E = 0.

    # R_CSM (R)
    if R_range[0] <= theta[3] <= R_range[-1]:
        prob_R = 1. / theta[3]
    else:
        prob_R = 0.

    # K_CSM (K)
    if K_range[0] <= theta[4] <= K_range[-1]:
        prob_K = 1. / theta[4]
    else:
        prob_K = 0.

    # Mixing (Mix)
    if Mix_range[0] <= theta[5] <= Mix_range[-1]:
        prob_Mix = 1. / theta[5]
    else:
        prob_Mix = 0.

    # Scaling factor (S)
    if S_range[0] <= theta[6] <= S_range[-1]:
        prob_S = 1. / theta[6]
    else:
        prob_S = 0.

    # T of explosion (T)
    if T_range[0] <= theta[7] <= T_range[-1]:
        prob_T = 1. / theta[7]
    else:
        prob_T = 0.

    # sum probabilities
    prob_total = np.log(prob_M * prob_Ni * prob_E * prob_R * prob_K * prob_Mix * prob_S * prob_T)
    return prob_total


def log_likelihood(theta, data):
    sampled = [Mzams_range, Ni_range, E_final_range, R_range, K_range, Mix_range]
    if (Mzams_range[0] <= theta[0] <= Mzams_range[-1]) & \
            (Ni_range[0] <= theta[1] <= Ni_range[-1]) & \
            (E_final_range[0] <= theta[2] <= E_final_range[-1]) & \
            (R_range[0] <= theta[3] <= R_range[-1]) & \
            (K_range[0] <= theta[4] <= K_range[-1]) & \
            (Mix_range[0] <= theta[5] <= Mix_range[-1]) & \
            (S_range[0] <= theta[6] <= S_range[-1]) & \
            (T_range[0] <= theta[7] <= T_range[-1]):
        log_likeli = 0
        data_x_allfilt = data['t_from_discovery'].unique() - 15 + theta[7]
        filters = list(data['filter'].unique())
        y_fit = interp.snec_interpolator(theta[0:6], sampled, data_x_allfilt, filters, pysynphot_models)
        if not isinstance(y_fit, str):
            # multiply whole graph by scaling factor
            y_fit = y_fit * theta[6]
            y_fit['time'] = data_x_allfilt
            for filt in filters:
                data_filt = data.loc[data['filter'] == filt]
                data_x_filt = data_filt['t_from_discovery'] - 15 + theta[7]
                data_y_filt = data_filt['abs_mag']
                data_dy_filt = data_filt['dmag']
                y_fit_filt = y_fit[filt].loc[y_fit['time'].isin(data_x_filt)]
                if not len(data_x_filt) == len(y_fit_filt):
                    print('stuck')
                log_likeli += -np.sum((data_y_filt - y_fit_filt) ** 2. / (2. * data_dy_filt ** 2.)
                                      + np.log(data_dy_filt))
                # print('log likelihood', log_likeli)
        else:
            log_likeli = - 10 ** 30  # just a very big number so it won't go past the edge values
            print('impossible SN')
    else:
        log_likeli = - 10 ** 30  # just a very big number so it won't go past the edge values
        print('not valid')
    # print('log likelihood', log_likeli)
    return log_likeli



def log_posterior(theta, data):
    ln_post = log_prior(theta) + log_likelihood(theta, data)
    return ln_post



def emcee_fit_params(data):
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=[data])

    Mzams_random = np.random.rand(n_walkers) * (Mzams_range[-1] - Mzams_range[0]) + Mzams_range[0]
    Ni_random = np.random.rand(n_walkers) * (Ni_range[-1] - Ni_range[0]) + Ni_range[0]
    E_random = np.random.rand(n_walkers) * (E_final_range[-1] - E_final_range[0]) + E_final_range[0]
    R_random = np.random.rand(n_walkers) * (R_range[-1] - R_range[0]) + R_range[0]
    K_random = np.random.rand(n_walkers) * (K_range[-1] - K_range[0]) + K_range[0]
    Mix_random = np.random.rand(n_walkers) * (Mix_range[-1] - Mix_range[0]) + Mix_range[0]
    S_random = np.random.rand(n_walkers) * (S_range[-1] - S_range[0]) + S_range[0]
    T_random = np.random.rand(n_walkers) * (T_range[-1] - T_range[0]) + T_range[0]

    initial_guesses = np.array([Mzams_random, Ni_random, E_random, R_random, K_random, Mix_random, S_random, T_random])
    initial_guesses = initial_guesses.T
    sampler.run_mcmc(initial_guesses, n_steps)

    return sampler


def chain_plots(sampler, **kwargs):
    chain = sampler.chain

    f_Mzams = plt.figure()
    plt.plot(chain[:, :, 0].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Mzams')
    plt.tight_layout()
    f_Mzams.savefig(os.path.join(res_dir, 'Mzams.png'))

    f_Ni = plt.figure()
    plt.plot(chain[:, :, 1].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Ni')
    plt.tight_layout()
    f_Ni.savefig(os.path.join(res_dir, 'Ni.png'))

    f_E = plt.figure()
    plt.plot(chain[:, :, 2].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('E')
    plt.tight_layout()
    f_E.savefig(os.path.join(res_dir, 'E.png'))

    f_R = plt.figure()
    plt.plot(chain[:, :, 3].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('R')
    plt.tight_layout()
    f_R.savefig(os.path.join(res_dir, 'R.png'))

    f_K = plt.figure()
    plt.plot(chain[:, :, 4].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('K')
    plt.tight_layout()
    f_K.savefig(os.path.join(res_dir, 'K.png'))

    f_Mix = plt.figure()
    plt.plot(chain[:, :, 5].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Mixing')
    plt.tight_layout()
    f_Mix.savefig(os.path.join(res_dir, 'Mix.png'))

    f_S = plt.figure()
    plt.plot(chain[:, :, 6].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Scaling')
    plt.tight_layout()
    f_S.savefig(os.path.join(res_dir, 'S.png'))

    f_T = plt.figure()
    plt.plot(chain[:, :, 7].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('T_exp')
    plt.tight_layout()
    f_T.savefig(os.path.join(res_dir, 'T.png'))


def get_each_walker_result(sampler, step):
    params = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
    dict = {}
    for i in range(len(params)):
        last_results = sampler.chain[:, step:, i]
        avg = np.average(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        # dict[params[i]+'_all_walkers'] = last_results
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg
    print(dict)

def get_param_results_dict(sampler, step):
    params = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
    dict = {}
    for i in range(len(params)):
        last_results = sampler.chain[:, step:, i]
        avg = np.average(last_results)
        if params[i] == 'T':
            avg -= 15
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        # dict[params[i]+'_all_walkers'] = last_results
        dict[params[i] + '_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg
    print(dict)

    with open(os.path.join(res_dir, 'final_results.csv'), 'w') as f:  # Just use 'w' mode in 3.x
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

def result_text_from_dict(sampler, step):
    param_dict = get_param_results_dict(sampler, step)
    res_text = SN_name + '\n'
    for param in ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']:
        res_text += param + ': ' + rounded_str(param_dict[param]) + r'$\pm$ [' +\
                        rounded_str(param_dict[param+'_lower']) + ',' +\
                        rounded_str(param_dict[param+'_upper']) + ']\n'
    return res_text






def plot_lightcurve_with_fit(SN_data, sampler, step):
    data_x = SN_data['t_from_discovery']
    filters = list(SN_data['filter'].unique())
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler.chain[i, step, :]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        data_x_moved_all = data_x - 15 + T
        sampled = [Mzams_range, Ni_range, E_final_range, R_range, K_range, Mix_range]
        y_fit = interp.snec_interpolator(requested[0:6], sampled, data_x_moved_all, filters, pysynphot_models)
        if not isinstance(y_fit, str):
            # multiply whole graph by scaling factor
            y_fit = y_fit * S
            for filt in filters:
                ax.plot(data_x_moved_all, y_fit[filt], color=colors[filt], alpha=0.3)
            log_likeli.append(log_likelihood(requested, SN_data))
    log_likeli = np.average(log_likeli)
    for filt in filters:
        data_filt = SN_data.loc[SN_data['filter'] == filt]
        data_x = data_filt['t_from_discovery']
        data_x_moved_filt = data_x + -15 + T
        data_y = data_filt['abs_mag']
        data_dy = data_filt['dmag']
        ax.errorbar(data_x_moved_filt, data_y, yerr=data_dy, marker='o', linestyle='None',
                    label=filt + ' ' + SN_name, color=colors[filt])
    results_text = result_text_from_dict(sampler, step)
    ax.text(0.6, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    ax.set_title('step '+str(step)
                 +'\nlog likelihood = ' + str(int(log_likeli)), fontsize=14)
    plt.tight_layout()
    ax.set_xlim(-2, 137)
    ax.set_ylim(14, 20)
    f_fit.savefig(os.path.join(res_dir, 'mag_lightcurve_fit' +str(step) +'.png'))
    return log_likeli



'''
# running code #
'''

sampler = emcee_fit_params(SN)

flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")



log_liklihood_mag = plot_lightcurve_with_fit(SN_w_early, sampler, n_steps-1)
# log_liklihood_mag = plot_lightcurve_with_fit(SN, sampler, 50)
# log_liklihood_mag = plot_lightcurve_with_fit(SN, sampler, 20)
log_liklihood_mag = plot_lightcurve_with_fit(SN_w_early, sampler, 0)

# to correct for T (time after explostion) actually being T+15
sampler.chain[:, :, 7] = - (sampler.chain[:, :, 7] - 15)
chain_plots(sampler)


labels = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
corner_range = [1., 1., 1., 1., 1., 1., 1., 1.]
f_corner = corner.corner(flat_sampler_no_burnin, labels=labels, range=corner_range)
# plt.tight_layout()
f_corner.savefig(os.path.join(res_dir, 'corner_plot.png'))

# MCMC_results = get_param_results_dict(sampler)

print(sampler.chain.shape)


