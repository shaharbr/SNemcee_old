import emcee
import numpy as np
from matplotlib import pyplot as plt
import snec_result_interpolator_mag as interp_mag
import snec_result_interpolator_lum as interp_lum
import snec_result_interpolator_veloc as interp_veloc
import pandas as pd
import datetime
from pathlib import Path
import os
import corner
import csv

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


# important to note: values can't be negative!
Mzams_range = [9.0, 11.0, 13.0, 15.0]
Ni_range = [0.02, 0.07, 0.12, 0.17]
E_final_range = [0.7, 1.0, 1.3, 1.9]
Mix_range = [2.0, 3.0]
R_range = [600, 3000]
K_range = [0.001, 120]
S_range = [0.5, 1.5]
T_range = [0, 30] # because can't have negative values, do 15 minus diff (so 0 is -15, and 30 is +15)
pysynphot_models = True

n_walkers = 80
n_steps = 51
n_params = 8
burn_in = 30


time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_combined')
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

# import SN2018hmx mag data
data_filepath = os.path.join('results', 'SN2018hmx_lightcurves')
data_mag = pd.read_csv(data_filepath, usecols=['dmag', 'filter', 'abs_mag', 't_from_discovery'])
if pysynphot_models:
    data_mag = data_mag.loc[data_mag['t_from_discovery'] < 140]
data_mag['abs_mag'] = data_mag['abs_mag'].abs()
filters = list(data_mag['filter'].unique())
filters = list(set(filters).intersection(['g', 'r', 'i']))
data_mag = data_mag.loc[data_mag['filter'].isin(filters)]
data_mag = data_mag.sort_values('t_from_discovery')
colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}

# import SN2018hmx bolometric lum data
data_lum = pd.read_csv(os.path.join('results', 'blackbody_results_18hmx_BVgri.csv'))
# convert watt to erg/s
data_lum['Lum'] = data_lum['Lum'] * 10**7
data_lum['dLum0'] = data_lum['dLum0'] * 10**7
data_lum['dLum1'] = data_lum['dLum1'] * 10**7


# replicate the last point of the photometry data x times to artificially increase weight of fitting to last point (10-folds)
times_to_amplify = 1
if times_to_amplify > 1:
    for data in [data_lum, data_mag]:
        last_row = data.loc[data['t_from_discovery'] > 350]
        last_row_repeats = pd.concat([last_row]*(times_to_amplify-1), ignore_index=True)
        data = pd.concat([data, last_row_repeats], ignore_index=True)

# import SN2018hmx photospheric velocities lum data
data_veloc = pd.read_csv(os.path.join('results', 'sN2018hmx_expansion_velocities.csv'),
                      usecols=['t_from_discovery', 'line', 'absorption_mean_velocity','absorption_std_velocity'])
data_veloc.rename({'absorption_mean_velocity':'veloc', 'absorption_std_velocity':'dveloc'}, axis='columns', inplace=True)
data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
# remove first point which seems like an artifact
data_veloc = data_veloc.loc[data_veloc['t_from_discovery'] > 20]



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
        filt = list(data['mag']['filter'].unique())
        # series of the times for each of the observations
        x_moved_mag = data['mag']['t_from_discovery'] - 15 + theta[7]
        x_moved_lum = data['lum']['t_from_discovery'] - 15 + theta[7]
        x_moved_veloc = data['veloc']['t_from_discovery'] - 15 + theta[7]
        # get the interpolated snec model for the parameter values in theta[0:6]
        y_fit_mag = interp_mag.snec_interpolator(theta[0:6], sampled, x_moved_mag, filt, pysynphot_models)
        y_fit_lum = interp_lum.snec_interpolator(theta[0:6], sampled, x_moved_lum)
        y_fit_veloc = interp_veloc.snec_interpolator(theta[0:6], sampled, x_moved_veloc)
        # multiply whole graph by scaling factor
        y_fit_lum = y_fit_lum * theta[6]
        y_fit_mag = y_fit_mag * theta[6]
        # add log likelihood for luminosity
        log_likelihood_lum = -np.sum((data['lum']['Lum'] - y_fit_lum) ** 2. /
                             (2. * data['lum']['dLum0'] ** 2.) + np.log(data['lum']['Lum']))
        log_likeli += log_likelihood_lum
        print('log likelihood mag', log_likelihood_lum)
        # add log likelihood for velocity
        log_likelihood_veloc = -np.sum((data['veloc']['veloc'] - y_fit_veloc) ** 2. /
                             (2. * data['veloc']['dveloc'] ** 2.) + np.log(data['veloc']['veloc']))
        log_likeli += log_likelihood_veloc
        print('log likelihood mag', log_likelihood_veloc)
        # add log likelihood for magnitudes
        # for the mag only, add the a time column so it can be filtered appropriately for each filter
        y_fit_mag['time'] = x_moved_mag
        log_likelihood_mag = 0
        # TODO yfilt needs to be a df with all filters and a column for time
        # data_filt = data.keys() #  TODO use this but ommit time column
        for filt in filters:
            data_filt = data['mag'].loc[data['mag']['filter'] == filt]
            data_filt_x = data_filt['t_from_discovery'] - 15 + theta[7]
            data_filt_y = data_filt['abs_mag']
            data_filt_dy = data_filt['dmag']
            y_fit_mag_filt = y_fit_mag[filt].loc[y_fit_mag['time'].isin(data_filt_x)]
            if not len(data_filt_x) == len(y_fit_mag_filt):
                print('stuck')
            log_likelihood_mag_filt = -np.sum((data_filt_y - y_fit_mag_filt) ** 2. / (2. * data_filt_dy ** 2.)
                                  + np.log(data_filt_y))
            log_likelihood_mag += log_likelihood_mag_filt
            log_likeli += log_likelihood_mag
            print('log likelihood mag' + str(filt), log_likelihood_mag_filt)
    else:
        log_likeli = 100000000000000
        print('not valid')
    print('likelihood_log total', log_likeli)
    print('likelihood total', np.exp(log_likeli))
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



def get_param_results_dict(sampler, step):
    params = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
    dict = {}
    for i in range(len(params)):
        last_results = sampler.chain[:, step:, i]
        avg = np.average(last_results)
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



def plot_lightcurve_with_fit(data, sampler, step):
    param_dict = get_param_results_dict(sampler, step)
    Mzams = param_dict['Mzams']
    Ni = param_dict['Ni']
    E = param_dict['E']
    R = param_dict['R']
    K = param_dict['K']
    Mix = param_dict['Mix']
    S = param_dict['S']
    T = param_dict['T']
    results_text = 'Mzams: ' + str(round(Mzams, 1)) + \
                   ' Ni: ' + str(round(Ni, 3)) + \
                   ' E: ' + str(round(E, 1)) + \
                   ' R: ' + str(int(R)) + \
                   ' K: ' + str(int(K)) + \
                   ' Mix: ' + str(round(Mix, 1)) + \
                   ' S: ' + str(round(S, 2)) + \
                   ' T: ' + str(round(T, 1))
    print(results_text)

    requested = [Mzams, Ni, E, R, K, Mix]
    sampled = [Mzams_range, Ni_range, E_final_range, R_range, K_range, Mix_range]
    # TODO here too - filters sould be got from function inputs

    log_likeli = 0
    filt = list(data['mag']['filter'].unique())
    # series of the times for each of the observations
    x_moved_mag = data['mag']['t_from_discovery'] + T
    x_moved_lum = data['lum']['t_from_discovery'] + T
    x_moved_veloc = data['veloc']['t_from_discovery'] + T
    # get the interpolated snec model for the parameter values in
    y_fit_mag = interp_mag.snec_interpolator(requested, sampled, x_moved_mag, filt, pysynphot_models)
    y_fit_lum = interp_lum.snec_interpolator(requested, sampled, x_moved_lum)
    y_fit_veloc = interp_veloc.snec_interpolator(requested, sampled, x_moved_veloc)
    # multiply whole graph by scaling factor
    y_fit_lum = y_fit_lum * S
    y_fit_mag = y_fit_mag * S
    # add log likelihood for luminosity
    log_likelihood_lum = -np.sum((data['lum']['Lum'] - y_fit_lum) ** 2. /
                                 (2. * data['lum']['dLum0'] ** 2.) + np.log(data['lum']['Lum']))
    log_likeli += log_likelihood_lum
    print('log likelihood mag', log_likelihood_lum)
    # add log likelihood for velocity
    log_likelihood_veloc = -np.sum((data['veloc']['veloc'] - y_fit_veloc) ** 2. /
                                   (2. * data['veloc']['dveloc'] ** 2.) + np.log(data['veloc']['veloc']))
    log_likeli += log_likelihood_veloc
    print('log likelihood mag', log_likelihood_veloc)
    # add log likelihood for magnitudes
    # for the mag only, add the a time column so it can be filtered appropriately for each filter
    y_fit_mag['time'] = x_moved_mag
    log_likelihood_mag_filt ={}
    log_likelihood_mag = 0
    f_mag, ax_mag = plt.subplots(figsize=(10, 8))
    for filt in filters:
        data_filt = data['mag'].loc[data['mag']['filter'] == filt]
        data_filt_x = data_filt['t_from_discovery'] + T
        data_filt_y = data_filt['abs_mag']
        data_filt_dy = data_filt['dmag']
        y_fit_mag_filt = y_fit_mag[filt].loc[y_fit_mag['time'].isin(data_filt_x)]
        if not len(data_filt_x) == len(y_fit_mag_filt):
            print('stuck')

        log_likelihood_mag_filt[filt] = -np.sum((data_filt_y - y_fit_mag_filt) ** 2. / (2. * data_filt_dy ** 2.)
                                          + np.log(data_filt_y))
        log_likelihood_mag += log_likelihood_mag_filt[filt]
        log_likeli += log_likelihood_mag
        print('log likelihood mag' + str(filt), log_likelihood_mag_filt)
        ax_mag.plot(data_filt_x, y_fit_mag_filt, color=colors[filt])
        ax_mag.errorbar(data_filt_x, data_filt_y, yerr=data_filt_dy,
                        marker='o', linestyle='None', label=filt + '_SN 2018hmx', color=colors[filt])
    print('likelihood_log total', log_likeli)
    print('likelihood total', np.exp(log_likeli))
    # mag plot
    ax_mag.set_title('step '+str(step) +
                     '\nTotal log likelihood = ' + str(int(log_likelihood_mag)) +
                     '\nFilter log likelihood = ' + str([(filt + ': ' +
                            str(int(log_likelihood_mag_filt[filt]))) for filt in filters]) +
                     '\nsnec model: ' + results_text, fontsize=14)
    plt.tight_layout()
    ax_mag.legend()
    f_mag.savefig(os.path.join(res_dir, 'lightcurve_fit_mag.png'))
    # velocity plot
    f_veloc, ax_veloc = plt.subplots(figsize=(10, 8))
    ax_veloc.errorbar(x_moved_veloc, data['veloc']['veloc'], yerr=data['veloc']['dveloc'],
                      marker='o', linestyle='None', label='SN 2018hmx')
    ax_veloc.plot(x_moved_veloc, y_fit_veloc, label='best fit:\n' + results_text)
    ax_veloc.legend()
    ax_veloc.set_title('step ' + str(step) +
                       '\nlog likelihood = ' + str(int(log_likelihood_veloc)), fontsize=14)
    plt.tight_layout()
    f_veloc.savefig(os.path.join(res_dir, 'lightcurve_fit_veloc.png'))
    # bolometric lum plot
    f_lum, ax_lum = plt.subplots(figsize=(10, 8))
    ax_lum.errorbar(x_moved_lum, data['lum']['Lum'],
                      yerr=[data['lum']['dLum0'], data['lum']['dLum1']],
                      marker='o', linestyle='None', label='SN 2018hmx')
    ax_lum.plot(x_moved_lum, y_fit_lum, label='best fit:\n' + results_text)
    ax_lum.legend()
    ax_lum.set_title('step ' + str(step) +
                       '\nlog likelihood = ' + str(int(log_likelihood_lum)), fontsize=14)
    plt.tight_layout()
    f_lum.savefig(os.path.join(res_dir, 'lightcurve_fit_lum.png'))
    # bolometric lum log plot
    f_lumlog, ax_lumlog = plt.subplots(figsize=(10, 8))
    ax_lumlog.errorbar(x_moved_lum, data['lum']['Lum'],
                    yerr=[data['lum']['dLum0'], data['lum']['dLum1']],
                    marker='o', linestyle='None', label='SN 2018hmx')
    ax_lumlog.plot(x_moved_lum, y_fit_lum, label='best fit:\n' + results_text)
    ax_lumlog.legend()
    ax_lumlog.set_title('step ' + str(step) +
                     '\nlog likelihood = ' + str(int(log_likelihood_lum)), fontsize=14)
    plt.tight_layout()
    ax_lumlog.set_yscale('log')
    f_lumlog.savefig(os.path.join(res_dir, 'lightcurve_fit_lum_log.png'))







SN_data_all = {'veloc': data_veloc, 'lum': data_lum, 'mag': data_mag}

sampler = emcee_fit_params(SN_data_all)

# to correct for T (time after explostion) actually being T+15
sampler.chain[:, :, 7] = sampler.chain[:, :, 7] - 15
chain_plots(sampler)

results_vec = plot_lightcurve_with_fit(SN_data_all, sampler, n_steps-1)


flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")


labels = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
corner_range = [1., 1., 1., 1., 1., 1., 1., 1.]
f_corner = corner.corner(flat_sampler_no_burnin, labels=labels, range=corner_range)
# plt.tight_layout()
f_corner.savefig(os.path.join(res_dir, 'corner_plot.png'))

# MCMC_results = get_param_results_dict(sampler)

print(sampler.chain.shape)

