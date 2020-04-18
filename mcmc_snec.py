import emcee
import numpy as np
from matplotlib import pyplot as plt
import snec_result_interpolator as interp
import pandas as pd
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


Mzams_range = [13.0, 16.0, 19.0]
Ni_range = [0.07, 0.10, 0.13, 0.16, 0.19]
E_final_range = [1.5, 1.8, 2.4]
Mix_range = [3.0]
R_range = [600, 2400]
K_range = [0.001, 30]

m_Solar = 1.989 * (10 ** 33)  # gram

sn18hmx = pd.read_csv(r'results\blackbody_results_18hmx.csv')

# convert watt to erg/s
sn18hmx['Lum'] = sn18hmx['Lum'] * 10**7
sn18hmx['dLum0'] = sn18hmx['dLum0'] * 10**7
sn18hmx['dLum1'] = sn18hmx['dLum1'] * 10**7


interp_days = np.arange(15, 170, 1)
interp_lum = np.interp(interp_days, sn18hmx['t_from_discovery'], sn18hmx['Lum'])
interp_dlum0 = np.interp(interp_days, sn18hmx['t_from_discovery'], sn18hmx['dLum0'])
interp_dlum1 = np.interp(interp_days, sn18hmx['t_from_discovery'], sn18hmx['dLum1'])

interp_sn18hmx = pd.DataFrame(
    {'t_from_discovery': interp_days,
     'Lum': interp_lum,
     'dLum0': interp_dlum0,
     'dLum1': interp_dlum1
    })



def log_prior(theta):
    # theta is a vector containing a specific set of parameters theta = [M, Ni, E, Mix, R, K]

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

    # Ni_boundary (Mix)
    if Mix_range[0] <= theta[3] <= Mix_range[-1]:
        prob_Mix = 1. / theta[3]
    else:
        prob_Mix = 0.

    # R_CSM (R)
    if R_range[0] <= theta[4] <= R_range[-1]:
        prob_R = 1. / theta[4]
    else:
        prob_R = 0.

    # K_CSM (K)
    if K_range[0] <= theta[5] <= K_range[-1]:
        prob_K = 1. / theta[5]
    else:
        prob_K = 0.

    # sum probabilities
    prob_total = np.log(prob_M * prob_Ni * prob_E * prob_Mix * prob_R * prob_K)
    return prob_total


def log_likelihood(theta, data_x, data_y, data_dy):
    sampled = [Mzams_range, Ni_range, E_final_range, Mix_range, R_range, K_range]
    if (Mzams_range[0] <= theta[0] <= Mzams_range[-1]) & (Ni_range[0] <= theta[1] <= Ni_range[-1])\
        & (E_final_range[0] <= theta[2] <= E_final_range[-1])& (Mix_range[0] <= theta[3] <= Mix_range[-1])\
        & (R_range[0] <= theta[4] <= R_range[-1])& (K_range[0] <= theta[5] <= K_range[-1]):
        y_fit = interp.snec_interpolator(theta, sampled)
        chi2 = (data_y - y_fit) ** 2. / (2. * data_dy ** 2.) + np.log(data_y)
    else:
        chi2 = 10000000000
    ln_like = -np.sum(chi2)
    print('chi_ln', ln_like)
    print('chi_exp', np.exp(ln_like))

    return ln_like



def log_posterior(theta, data_x, data_y, data_dy):
    ln_post = log_prior(theta) + log_likelihood(theta, data_x, data_y, data_dy)
    return ln_post



def emcee_fit_params(data_time, data_lum, data_dlum):
    n_walkers = 30
    n_steps = 100
    n_params = 6
    args = [data_time, data_lum, data_dlum]
    # TODO only using dLum0, not dLum1?
    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=args)

    Mzams_random = np.random.rand(n_walkers) * (Mzams_range[-1] - Mzams_range[0]) + Mzams_range[0]
    Ni_random = np.random.rand(n_walkers) * (Ni_range[-1] - Ni_range[0]) + Ni_range[0]
    E_random = np.random.rand(n_walkers) * (E_final_range[-1] - E_final_range[0]) + E_final_range[0]
    Mix_random = np.random.rand(n_walkers) * (Mix_range[-1] - Mix_range[0]) + Mix_range[0]
    R_random = np.random.rand(n_walkers) * (R_range[-1] - R_range[0]) + R_range[0]
    K_random = np.random.rand(n_walkers) * (K_range[-1] - K_range[0]) + K_range[0]

    initial_guesses = np.array([Mzams_random, Ni_random, E_random, Mix_random, R_random, K_random])
    initial_guesses = initial_guesses.T

    sampler.run_mcmc(initial_guesses, n_steps)

    return sampler


# def get_LCO_V_df(SN_dict):
#     LCO_lighcurve = SN_dict['lightcurve']['Las Cumbres']['df']
#     V_LCO_lightcurve = LCO_lighcurve.loc[LCO_lighcurve['filter'] == 'V']
#     return V_LCO_lightcurve


def SN_lightcurve_params(SN_data):
    data_time = SN_data['t_from_discovery']
    data_Lum = SN_data['Lum']
    data_dLum = SN_data['dLum0']
    sampler = emcee_fit_params(data_time, data_Lum, data_dLum)
    return sampler


def chain_plots(chain, **kwargs):
    plt.figure()
    plt.plot(chain[:, :, 0].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Mzams')

    plt.figure()
    plt.plot(chain[:, :, 1].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Ni')

    plt.figure()
    plt.plot(chain[:, :, 2].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('E')

    plt.figure()
    plt.plot(chain[:, :, 3].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Mix')

    plt.figure()
    plt.plot(chain[:, :, 4].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('R')

    plt.figure()
    plt.plot(chain[:, :, 5].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('K')

# [walkers, step, dim]


def get_param_results_dict(sampler):
    params = ['Mzams', 'Ni', 'E', 'Mix', 'R', 'K']
    dict = {}
    for i in range(len(params)):
        last_results = sampler.chain[:, -1:, i]
        avg = np.average(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        dict[params[i]+'_all_walkers'] = last_results
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg
    print(dict)
    return dict



def calc_chi_square_sampled(data, x_fit, y_fit):
    sampling = np.arange(15, 170, 5)
    data_sampled = np.interp(sampling, data['t_from_discovery'], data['Lum'])
    data_err_sampled = np.interp(sampling, data['t_from_discovery'], data['dLum0'])
    model_sampled = np.interp(sampling, x_fit, y_fit)
    chisq = np.sum(((data_sampled - model_sampled) /
                    data_err_sampled) ** 2)
    chisq_reduced = chisq / (len(sampling) - 1)
    plt.figure()
    plt.plot(sampling, data_sampled, marker='o')
    plt.plot(sampling, model_sampled, marker='o')
    return chisq_reduced


def plot_lightcurve_with_fit(SN_data, sampler):
    param_dict = get_param_results_dict(sampler)
    Mzams = param_dict['Mzams']
    Ni = param_dict['Ni']
    E = param_dict['E']
    Mix = param_dict['Mix']
    R = param_dict['R']
    K = param_dict['K']
    results_text = 'Mzams: '+str(round(Mzams, 1))+' Ni: '+str(round(Ni, 3))+' E: '+str(round(E, 1))+' Mix: '+str(Mix)+' R: '+str(int(R))+' K: '+str(int(K))
    print(results_text)
    x_fit = interp_days
    data_x = SN_data['t_from_discovery']
    data_y = SN_data['Lum']
    dy0 = SN_data['dLum0']
    dy1 = SN_data['dLum1']

    requested = [Mzams, Ni, E, Mix, R, K]
    sampled = [Mzams_range, Ni_range, E_final_range, Mix_range, R_range, K_range]
    y_fit = interp.snec_interpolator(requested, sampled)

    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    ax.errorbar(data_x, data_y, yerr=[dy0, dy1], marker='o', linestyle='None', label='SN 2018hmx')
    ax.plot(x_fit, y_fit, label='best fit:\n'+results_text)
    ax.legend()
    chisq = calc_chi_square_sampled(SN_data, x_fit, y_fit)
    ax.set_title('chi_sq_red = ' + str(int(chisq)), fontsize=14)





sampler = SN_lightcurve_params(interp_sn18hmx)
chain_plots(sampler.chain)
results_vec = plot_lightcurve_with_fit(sn18hmx, sampler)


MCMC_results = get_param_results_dict(sampler)
# Ni_results['Ni_87A'] = np.average([Ni_mass_by_slope(i, line, SN87A_line) for i in [150, 300]])
# print(Ni_results)
# param_results = pd.DataFrame(Ni_results, index=[0])
# param_results.to_csv(r'results\Ni_results_' + SN + '_BVgri.csv')

print(sampler.chain.shape)



fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)


# Mzams = 10
# Ni = 0.1
# E = 1.1
# Mix = 3
# R = 500
# K = 20
# results_text = 'Mzams: '+str(round(Mzams, 1))+' Ni: '+str(round(Ni, 3))+' E: '+str(round(E, 1))+' Mix: '+str(Mix)+' R: '+str(int(R))+' K: '+str(int(K))
# print(results_text)
# x_fit = interp_days
# data_x = sn18hmx['t_from_discovery']
# data_y = sn18hmx['Lum']
# dy0 = sn18hmx['dLum0']
# dy1 = sn18hmx['dLum1']

# requested = [Mzams, Ni, E, Mix, R, K]
# sampled = [Mzams_range, Ni_range, E_final_range, Mix_range, R_range, K_range]
# y_fit = interp.snec_interpolator(requested, sampled)

# ax.errorbar(data_x, data_y, yerr=[dy0, dy1], marker='o', linestyle='None', label='SN 2018hmx')
# ax.plot(x_fit, y_fit, label='best fit:\n'+results_text)
# ax.legend()
# chisq = calc_chi_square_sampled(sn18hmx, x_fit, y_fit)
# ax.set_title('chi_sq_red = ' + str(int(chisq)), fontsize=14)


# E = 2.0
# results_text = 'Mzams: '+str(round(Mzams, 1))+' Ni: '+str(round(Ni, 3))+' E: '+str(round(E, 1))+' Mix: '+str(Mix)+' R: '+str(int(R))+' K: '+str(int(K))
# print(results_text)

# requested = [Mzams, Ni, E, Mix, R, K]
# sampled = [Mzams_range, Ni_range, E_final_range, Mix_range, R_range, K_range]
# y_fit = interp.snec_interpolator(requested, sampled)
#
# ax.errorbar(data_x, data_y, yerr=[dy0, dy1], marker='o', linestyle='None', label='SN 2018hmx')
# ax.plot(x_fit, y_fit, label='best fit:\n'+results_text)
# ax.legend()
# chisq = calc_chi_square_sampled(sn18hmx, x_fit, y_fit)
# ax.set_title('chi_sq_red = ' + str(int(chisq)), fontsize=14)


plt.show()