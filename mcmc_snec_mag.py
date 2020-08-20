import emcee
import numpy as np
from matplotlib import pyplot as plt
import snec_result_interpolator_mag as interp
import pandas as pd
import datetime
from pathlib import Path
import os
import corner




'''
parts of this code are based on code by Griffin Hosseinzadeh
'''
# TODO run MCMC
# TODO run mmc with more burn in + 350, wlaker 50. remove constant parameters
#

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



data_filepath = os.path.join('results', 'SN2018hmx_lightcurves')
SN = pd.read_csv(data_filepath, usecols=['dmag', 'filter', 'abs_mag', 't_from_discovery'])
print(SN)
SN['abs_mag'] = SN['abs_mag'].abs()

filters = list(SN['filter'].unique())
filters = list(set(filters).intersection(['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']))
print(filters)
SN = SN.loc[SN['filter'].isin(filters)]
SN = SN.sort_values('t_from_discovery')
colors = {'u': 'purple', 'g': 'teal', 'r': 'red', 'i': 'maroon', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'green', 'R': 'red', 'I': 'maroon'}



Mzams_range = [13.0, 16.0, 19.0, 21.0]
Ni_range = [0.10, 0.13, 0.16, 0.19]
E_final_range = [1.5, 1.8, 2.4]
R_range = [600, 1800, 2400, 3000]
K_range = [0.001, 30, 90]


n_walkers = 10
n_steps = 20
n_params = 5

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
Path(os.path.join('mcmc_results', str(time_now))).mkdir(parents=True, exist_ok=True)

run_param_df = pd.DataFrame({'parameter_ranges': {'Mzams_range': Mzams_range, 'Ni_range': Ni_range,
                                                  'E_final_range': E_final_range,
                                                  'R_range': R_range, 'K_range': K_range},
                             'MCMC_parameters': {'n_walkers': n_walkers, 'n_steps': n_steps, 'n_params': n_params},
                             'time': time_now})
run_param_df.to_csv(os.path.join('mcmc_results', str(time_now), 'run_parameters.csv'))


m_Solar = 1.989 * (10 ** 33)  # gram





def log_prior(theta):
    # theta is a vector containing a specific set of parameters theta = [M, Ni, E, R, K]

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

    # sum probabilities
    prob_total = np.log(prob_M * prob_Ni * prob_E * prob_R * prob_K)
    return prob_total


def log_likelihood(theta, data):
    sampled = [Mzams_range, Ni_range, E_final_range, R_range, K_range]
    if (Mzams_range[0] <= theta[0] <= Mzams_range[-1]) & (Ni_range[0] <= theta[1] <= Ni_range[-1])\
        & (E_final_range[0] <= theta[2] <= E_final_range[-1])\
        & (R_range[0] <= theta[3] <= R_range[-1])& (K_range[0] <= theta[4] <= K_range[-1]):
        chisq_ln_like = 0
        data_x_allfilt = data['t_from_discovery']
        # TODO yfilt needs to be a df with all filters and a column for time
        data_filt = data.keys() #  TODO use this but ommit time column
        # TODO filters here is taken from outside function - should get it in inputs to the function (which should get it from upstream funcitons)
        y_fit = interp.snec_interpolator(theta, sampled, data_x_allfilt, filters)
        y_fit['time'] = data_x_allfilt
        for filt in filters:
            print(filt)
            data_filt = data.loc[data['filter'] == filt]
            data_x = data_filt['t_from_discovery']
            data_y = data_filt['abs_mag']
            data_dy = data_filt['dmag']
            # y_fit = y_fit.reset_index(drop=True)
            # y_fit = y_fit.rename(columns={'time': 't_from_discovery'})
            y_fit_filt = y_fit[filt].loc[y_fit['time'].isin(data_x)]
            # y_fit_filt = y_fit.merge(data_x, on='t_from_discovery', how='right')
            # print('merge', data.merge(data_x, on='t_from_discovery', how='right'))
            # y_fit_filt = y_fit.loc[y_fit['time'].astype(str).contains('|'.join(data_x.astype(str)))]
            if not len(data_x) == len(y_fit_filt):
                print('stuck')
                exit()
            # TODO is the chisq here correct?
            print('a0', data_y)
            print('a1', data_y - y_fit_filt)
            print('a2', (data_y - y_fit_filt) ** 2.)
            print('a3', (2. * data_dy ** 2.))
            print('a4', np.log(data_y))
            chisq_ln_like += -np.sum((data_y - y_fit_filt) ** 2. / (2. * data_dy ** 2.) + np.log(data_y))
            print('chi2', chisq_ln_like)
    else:
        chisq_ln_like = 100000000000000

    print('chi_ln', chisq_ln_like)
    print('chi_exp', np.exp(chisq_ln_like))
    return chisq_ln_like



def log_posterior(theta, data):
    ln_post = log_prior(theta) + log_likelihood(theta, data)
    return ln_post



def emcee_fit_params(data):
    # TODO only using dLum0, not dLum1?
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=[data])

    Mzams_random = np.random.rand(n_walkers) * (Mzams_range[-1] - Mzams_range[0]) + Mzams_range[0]
    Ni_random = np.random.rand(n_walkers) * (Ni_range[-1] - Ni_range[0]) + Ni_range[0]
    E_random = np.random.rand(n_walkers) * (E_final_range[-1] - E_final_range[0]) + E_final_range[0]
    R_random = np.random.rand(n_walkers) * (R_range[-1] - R_range[0]) + R_range[0]
    K_random = np.random.rand(n_walkers) * (K_range[-1] - K_range[0]) + K_range[0]

    initial_guesses = np.array([Mzams_random, Ni_random, E_random, R_random, K_random])
    initial_guesses = initial_guesses.T

    sampler.run_mcmc(initial_guesses, n_steps)

    return sampler


# def get_LCO_V_df(SN_dict):
#     LCO_lighcurve = SN_dict['lightcurve']['Las Cumbres']['df']
#     V_LCO_lightcurve = LCO_lighcurve.loc[LCO_lighcurve['filter'] == 'V']
#     return V_LCO_lightcurve


def SN_lightcurve_params(SN_data):
    # TODO remove, redundant
    # data_time = SN_data['t_from_discovery']
    # data_Lum = SN_data['Lum']
    # data_dLum = SN_data['dLum0']
    sampler = emcee_fit_params(SN_data)
    return sampler


def chain_plots(sampler, **kwargs):
    chain = sampler.chain

    f_Mzams = plt.figure()
    plt.plot(chain[:, :, 0].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Mzams')
    plt.tight_layout()
    f_Mzams.savefig(os.path.join('mcmc_results', str(time_now), 'Mzams.png'))

    f_Ni = plt.figure()
    plt.plot(chain[:, :, 1].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Ni')
    plt.tight_layout()
    f_Ni.savefig(os.path.join('mcmc_results', str(time_now), 'Ni.png'))

    f_E = plt.figure()
    plt.plot(chain[:, :, 2].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('E')
    plt.tight_layout()
    f_E.savefig(os.path.join('mcmc_results', str(time_now), 'E.png'))

    f_R = plt.figure()
    plt.plot(chain[:, :, 3].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('R')
    plt.tight_layout()
    f_R.savefig(os.path.join('mcmc_results', str(time_now), 'R.png'))

    f_K = plt.figure()
    plt.plot(chain[:, :, 4].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('K')
    plt.tight_layout()
    f_K.savefig(os.path.join('mcmc_results', str(time_now), 'K.png'))

# [walkers, step, dim]

import csv

def get_param_results_dict(sampler, step):
    params = ['Mzams', 'Ni', 'E', 'R', 'K']
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

    with open(os.path.join('mcmc_results', str(time_now), 'final_results.csv'), 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)


    # df = pd.DataFrame.from_dict(data=dict, orient='index', index=)
    # df.to_csv(os.path.join('mcmc_results', str(time_now), 'final_results.csv'))
    return dict



def plot_chi2(SN_filt_data, y_fit, filter):
    f_chi = plt.figure()
    plt.plot(SN_filt_data['t_from_discovery'], SN_filt_data['abs_mag'], marker='o')
    plt.plot(SN_filt_data['t_from_discovery'], y_fit, marker='o')
    plt.tight_layout()
    f_chi.savefig(os.path.join('mcmc_results', str(time_now), 'chi_square_sampling_'+filter+'.png'))


def plot_lightcurve_with_fit(SN_data, sampler, step):
    param_dict = get_param_results_dict(sampler, step)
    Mzams = param_dict['Mzams']
    Ni = param_dict['Ni']
    E = param_dict['E']
    R = param_dict['R']
    K = param_dict['K']
    results_text = 'Mzams: '+str(round(Mzams, 1))+' Ni: '+str(round(Ni, 3))+' E: '+str(round(E, 1))+' R: '+str(int(R))+' K: '+str(int(K))
    print(results_text)
    data_x = SN_data['t_from_discovery']
    data_y = SN_data['abs_mag']
    dy = SN_data['dmag']

    requested = [Mzams, Ni, E, R, K]
    sampled = [Mzams_range, Ni_range, E_final_range, R_range, K_range]
    # TODO here too - filters sould be got from function inputs
    data_x_allfilt = SN_data['t_from_discovery']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    # TODO here too filters
    chi2 = 0
    for filt in filters:
        data_filt = SN_data.loc[SN_data['filter'] == filt]
        data_x = data_filt['t_from_discovery']
        data_y = data_filt['abs_mag']
        data_dy = data_filt['dmag']
        y_fit_filt = interp.snec_interpolator(requested, sampled, data_x, [filt])[filt]
        chi2 += np.sum(((data_y - y_fit_filt) /data_dy)**2)
        ax.plot(data_x, y_fit_filt, label=filt+' best fit:\n' + results_text, color=colors[filt])
        ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', label=filt + '_SN 2018hmx', color=colors[filt])
        # plot_chi2(data_filt, y_fit_filt, filt)
    chi2_reduced = chi2 / (len(SN_data['t_from_discovery']) - 1)
    ax.set_title('step '+str(step)+'\nchi_sq_red = ' + str(int(chi2_reduced)), fontsize=14)
    plt.tight_layout()
    ax.legend()
    f_fit.savefig(os.path.join('mcmc_results', str(time_now), 'lightcurve_fit.png'))





sampler = SN_lightcurve_params(SN)
chain_plots(sampler)
results_vec = plot_lightcurve_with_fit(SN, sampler, 0)
results_vec = plot_lightcurve_with_fit(SN, sampler, 19)


flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join('mcmc_results', str(time_now), 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=0, flat=True)
np.savetxt(os.path.join('mcmc_results', str(time_now), 'flat_sampler_without100burnin.csv'), flat_sampler_no_burnin, delimiter=",")


labels = ['Mzams', 'Ni', 'E', 'R', 'K']
corner_range = [1., 1., 1., 1., 1.,]
f_corner = corner.corner(flat_sampler_no_burnin, labels=labels, range=corner_range)
# plt.tight_layout()
f_corner.savefig(os.path.join('mcmc_results', str(time_now), 'corner_plot.png'))

# MCMC_results = get_param_results_dict(sampler)

print(sampler.chain.shape)

plt.show()

