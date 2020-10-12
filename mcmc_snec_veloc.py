import emcee
import numpy as np
from matplotlib import pyplot as plt
import snec_result_interpolator_veloc as interp
import pandas as pd
import datetime
from pathlib import Path
import os
import corner
import csv
# TODO add mixing, scale and time from explosion parameters
# TODO make two codes, one with CSM and one without? (to compare to other mcmc paper)

'''
parts of this code are based on code by Griffin Hosseinzadeh
'''
# TODO run MCMC
# TODO run mmc with more burn in + 350, wlaker 50. remove constant parameters
# TODO document code
# TODO do the same for the mag code
#

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



# important to note: values can't be negative!
Mzams_range = [12.0, 15.0]
Ni_range = [0.02, 0.07, 0.12, 0.17]
E_final_range = [1.2, 1.7, 2.2, 2.7]
Mix_range = [1.0, 3.0, 6.0]
R_range = [600, 1400, 2200, 3000]
K_range = [0.001, 50, 100, 150]
S_range = [0.5, 2.0]
T_range = [0, 30] # because can't have negative values, do 15 minus diff (so 0 is -15, and 30 is +15)



n_walkers = 16
n_steps = 30
n_params = 8
burn_in = 0

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
Path(os.path.join('mcmc_results', str(time_now)+'_veloc')).mkdir(parents=True, exist_ok=True)

run_param_df = pd.DataFrame.from_dict({'Mzams_range': str(Mzams_range), 'Ni_range': str(Ni_range),
                             'E_final_range': str(E_final_range), 'R_range': str(R_range), 'K_range': str(K_range),
                             'Mix_range': str(Mix_range), 'Scaling_range': str(S_range), 'T_exp_min15_range': str(T_range),
                             'n_walkers': n_walkers, 'n_steps': n_steps, 'n_params': n_params,
                             'burn_in': burn_in,
                             'time': time_now}, orient='index')
run_param_df.to_csv(os.path.join('mcmc_results', str(time_now), 'run_parameters.csv'))


m_Solar = 1.989 * (10 ** 33)  # gram

sn18hmx = pd.read_csv(os.path.join('results', 'sN2018hmx_expansion_velocities.csv'),
                      usecols=['t_from_discovery', 'line', 'absorption_mean_velocity','absorption_std_velocity'])
sn18hmx.rename({'absorption_mean_velocity':'veloc', 'absorption_std_velocity':'dveloc'}, axis='columns', inplace=True)
sn18hmx = sn18hmx.loc[sn18hmx['line'] == 'FeII 5169']
sn18hmx.sort_values(by=['t_from_discovery'], inplace=True)
# remove first point which seems like an artifact
sn18hmx = sn18hmx.loc[sn18hmx['t_from_discovery'] > 20]
# TODO this is a patch until i calc the real dy on late velocities, remember to fix!
sn18hmx['dveloc'][35] = 80

# plt.figure()
# plt.plot(sn18hmx['t_from_discovery'], sn18hmx['veloc'], marker='o')


# remove first timepoint
# sn18hmx = sn18hmx.loc[sn18hmx['t_from_discovery'] > 20]

# replicate the last point x times to artificially increase weight of fitting to last point (10-folds)
times_to_amplify = 1
if times_to_amplify > 1:
    last_row = sn18hmx.loc[sn18hmx['t_from_discovery'] > 350]
    last_row_repeats = pd.concat([last_row]*(times_to_amplify-1), ignore_index=True)
    sn18hmx = pd.concat([sn18hmx, last_row_repeats], ignore_index=True)


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


def log_likelihood(theta, data_x, data_y, data_dy):
    sampled = [Mzams_range, Ni_range, E_final_range, R_range, K_range, Mix_range]
    if (Mzams_range[0] <= theta[0] <= Mzams_range[-1]) &\
            (Ni_range[0] <= theta[1] <= Ni_range[-1]) &\
            (E_final_range[0] <= theta[2] <= E_final_range[-1]) &\
            (R_range[0] <= theta[3] <= R_range[-1]) &\
            (K_range[0] <= theta[4] <= K_range[-1]) & \
            (Mix_range[0] <= theta[5] <= Mix_range[-1]) & \
            (S_range[0] <= theta[6] <= S_range[-1]) & \
            (T_range[0] <= theta[7] <= T_range[-1]):
        data_x_moved = data_x - 15 + theta[7]
        y_fit = interp.snec_interpolator(theta[0:6], sampled, data_x_moved)
        # multiply whole graph by scaling factor
        y_fit = y_fit * theta[6]

        chi2 = (data_y - y_fit) ** 2. / (2. * data_dy ** 2.) + np.log(data_y)
    else:
        chi2 = 10000000000
    ln_like = -np.sum(chi2)
    print('chi_ln', ln_like)
    return ln_like


def log_posterior(theta, data_x, data_y, data_dy):
    ln_post = log_prior(theta) + log_likelihood(theta, data_x, data_y, data_dy)
    return ln_post



def emcee_fit_params(data_time, data_veloc, data_dveloc):
    args = [data_time, data_veloc, data_dveloc]
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=args)

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


# def get_LCO_V_df(SN_dict):
#     LCO_lighcurve = SN_dict['lightcurve']['Las Cumbres']['df']
#     V_LCO_lightcurve = LCO_lighcurve.loc[LCO_lighcurve['filter'] == 'V']
#     return V_LCO_lightcurve


def SN_lightcurve_params(SN_data):
    data_time = SN_data['t_from_discovery']
    data_veloc = SN_data['veloc']
    data_dveloc = SN_data['dveloc']
    sampler = emcee_fit_params(data_time, data_veloc, data_dveloc)
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

    f_Mix = plt.figure()
    plt.plot(chain[:, :, 5].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Mixing')
    plt.tight_layout()
    f_Mix.savefig(os.path.join('mcmc_results', str(time_now), 'Mix.png'))

    f_S = plt.figure()
    plt.plot(chain[:, :, 6].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Scaling')
    plt.tight_layout()
    f_S.savefig(os.path.join('mcmc_results', str(time_now), 'S.png'))

    f_T = plt.figure()
    plt.plot(chain[:, :, 7].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('T_exp')
    plt.tight_layout()
    f_T.savefig(os.path.join('mcmc_results', str(time_now), 'T.png'))



def get_param_results_dict(sampler, step):
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

    with open(os.path.join('mcmc_results', str(time_now), 'final_results.csv'), 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)

    # df = pd.DataFrame.from_dict(data=dict, orient='index', index=)
    # df.to_csv(os.path.join('mcmc_results', str(time_now), 'final_results.csv'))
    return dict



def calc_chi_square_sampled(data, y_fit):
    # TODO is the chisq formula here correct? should be np.sum(((data_y - y_fit_filt)**2 /data_y)) ?
    chisq = np.sum(((data['veloc'] - y_fit) /
                    data['dveloc']) ** 2)
    chisq_reduced = chisq / (len(data['t_from_discovery']) - 1)
    f_chi = plt.figure()
    plt.plot(data['t_from_discovery'], data['veloc'], marker='o')
    plt.plot(data['t_from_discovery'], y_fit, marker='o')
    plt.tight_layout()
    f_chi.savefig(os.path.join('mcmc_results', str(time_now), 'chi_square_sampling.png'))
    return chisq_reduced


def plot_lightcurve_with_fit(SN_data, sampler, step):
    param_dict = get_param_results_dict(sampler, step)
    Mzams = param_dict['Mzams']
    Ni = param_dict['Ni']
    E = param_dict['E']
    R = param_dict['R']
    K = param_dict['K']
    Mix = param_dict['Mix']
    S = param_dict['S']
    T = param_dict['T']
    results_text = 'Mzams: '+str(round(Mzams, 1))+\
                   ' Ni: '+str(round(Ni, 3))+\
                   ' E: '+str(round(E, 1))+\
                   ' R: '+str(int(R))+\
                   ' K: '+str(int(K))+\
                   ' Mix: '+str(round(Mix, 1))+\
                   ' S: '+str(round(S, 2))+\
                   ' T: '+str(round(T, 1))
    print(results_text)
    data_x = SN_data['t_from_discovery']
    data_x_moved =  data_x + T
    data_y = SN_data['veloc']
    dy = SN_data['dveloc']

    requested = [Mzams, Ni, E, R, K, Mix]
    sampled = [Mzams_range, Ni_range, E_final_range, R_range, K_range, Mix_range]
    y_fit = interp.snec_interpolator(requested, sampled, data_x_moved)
    # multiply whole graph by scaling factor
    y_fit = y_fit * S

    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.errorbar(data_x_moved, data_y, yerr=dy, marker='o', linestyle='None', label='SN 2018hmx')
    ax.plot(data_x_moved, y_fit, label='best fit:\n'+results_text)
    ax.legend()
    chisq = calc_chi_square_sampled(SN_data, y_fit)
    ax.set_title('step '+str(step)+'\nchi_sq_red = ' + str(int(chisq)), fontsize=14)
    plt.tight_layout()
    f_fit.savefig(os.path.join('mcmc_results', str(time_now), 'lightcurve_fit.png'))



'''
# running code #
'''

sampler = SN_lightcurve_params(sn18hmx)
# to correct for T (time after explostion) actually being T+15
sampler.chain[:, :, 7] = sampler.chain[:, :, 7] - 15
chain_plots(sampler)
# results_vec = plot_lightcurve_with_fit(sn18hmx, sampler, 0)
# results_vec = plot_lightcurve_with_fit(sn18hmx, sampler, 3)
results_vec = plot_lightcurve_with_fit(sn18hmx, sampler, n_steps-1)


flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join('mcmc_results', str(time_now), 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
np.savetxt(os.path.join('mcmc_results', str(time_now), 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")


labels = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
corner_range = [1., 1., 1., 1., 1., 1., 1., 1.]
f_corner = corner.corner(flat_sampler_no_burnin, labels=labels, range=corner_range)
# plt.tight_layout()
f_corner.savefig(os.path.join('mcmc_results', str(time_now), 'corner_plot.png'))

# MCMC_results = get_param_results_dict(sampler, 0)
# Ni_results['Ni_87A'] = np.average([Ni_mass_by_slope(i, line, SN87A_line) for i in [150, 300]])
# print(Ni_results)
# param_results = pd.DataFrame(Ni_results, index=[0])
# param_results.to_csv(r'results\Ni_results_' + SN + '_BVgri.csv')

print(sampler.chain.shape)

plt.show()
