import emcee
import numpy as np
from matplotlib import pyplot as plt

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

Mzams_range = [10, 20]
Ni_mass_range = [0.1, 0.2]
E_final_range = [1, 2]
Ni_boundary_range = [3, 6]
R_CSM_range = [500, 10000]
K_CSM_range = [0.001, 100]

m_Solar = 1.989 * (10 ** 33)  # gram

def log_prior(theta):
    # theta is a vector containing a specific set of parameters theta = [M, Ni, E, Mix, R, K]

    # Mzams (M)
    if min(Mzams_range) < theta[0] < max(Mzams_range):
        prob_M = 1. / theta[0]
    else:
        prob_M = 0.

    # Ni mass (Ni)
    if min(Ni_mass_range) < theta[1] < max(Ni_mass_range):
        prob_Ni = 1. / theta[1]
    else:
        prob_Ni = 0.

    # E_final (E)
    if min(E_final_range) < theta[2] < max(E_final_range):
        prob_E = 1. / theta[2]
    else:
        prob_E = 0.

    # Ni_boundary (Mix)
    if min(Ni_boundary_range) < theta[3] < max(Ni_boundary_range):
        prob_Mix = 1. / theta[3]
    else:
        prob_Mix = 0.

    # R_CSM (R)
    if min(R_CSM_range) < theta[4] < max(R_CSM_range):
        prob_R = 1. / theta[4]
    else:
        prob_R = 0.

    # K_CSM (K)
    if min(K_CSM_range) < theta[5] < max(K_CSM_range):
        prob_K = 1. / theta[5]
    else:
        prob_K = 0.

    prob_total = np.log(prob_M * prob_Ni * prob_E * prob_Mix * prob_R * prob_K)
    return prob_total

def log_likelihood(theta, data_x, data_y, data_dy):
    [M, Ni, E, Mix, R, K] = theta

    # TODO:
    '''
    for each parameter:
    find closest two graphs 
    find closest two graphs 
    '''

    y_fit = M + p0 * data_x - a0 / (1 + np.exp((data_x - tPT) / w0))
    chi2 = (data_y - y_fit) ** 2. / (2. * data_dy ** 2.) + np.log(data_y)
    ln_like = -np.sum(chi2)
    return ln_like



def log_posterior(theta, data_x, data_y, data_dy):
    ln_post = log_prior(theta) + log_likelihood(theta, data_x, data_y, data_dy)
    return ln_post



def emcee_fit_params(v_time_vec, v_mag_vec, v_dmag_vec):
    n_walkers = 20
    n_steps = 300
    n_params = 4
    args = [v_time_vec, v_mag_vec, v_dmag_vec]

    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=args)

    tPT_min = 50
    tPT_max = 150
    a0_min = 0
    a0_max = 6
    w0_min = 0
    w0_max = 100
    m0_min = 10
    m0_max = 25

    m0_random = np.random.rand(n_walkers) * (m0_max - m0_min) + m0_min
    a0_random = np.random.rand(n_walkers) * (a0_max - a0_min) + a0_min
    tPT_random = np.random.rand(n_walkers) * (tPT_max- tPT_min) + tPT_min
    w0_random = np.random.rand(n_walkers) * (w0_max - w0_min) + w0_min
    initial_guesses = np.array([m0_random, a0_random, tPT_random, w0_random]).T

    sampler.run_mcmc(initial_guesses, n_steps)

    return sampler


def get_LCO_V_df(SN_dict):
    LCO_lighcurve = SN_dict['lightcurve']['Las Cumbres']['df']
    V_LCO_lightcurve = LCO_lighcurve.loc[LCO_lighcurve['filter'] == 'V']
    return V_LCO_lightcurve


def SN_lightcurve_params(SN_dict):
    df = get_LCO_V_df(SN_dict)
    # replicate the dots at the end of the light curve
    v_time_vec = df['t_from_discovery']
    v_mag_vec = df['mag']
    v_dmag_vec = df['dmag']
    sampler = emcee_fit_params(v_time_vec, v_mag_vec, v_dmag_vec)
    return sampler


def chain_plots(chain, **kwargs):
    plt.figure()
    plt.plot(chain[:, :, 0].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('m0')

    plt.figure()
    plt.plot(chain[:, :, 1].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('a0')

    plt.figure()
    plt.plot(chain[:, :, 2].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('tPT')

    plt.figure()
    plt.plot(chain[:, :, 3].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('w0')

# [walkers, step, dim]


def get_param_results_dict(sampler):
    params = ['m0', 'a0', 'tPT', 'w0']
    dict = {}
    for i in range(len(params)):
        last_results = sampler.chain[:, -100:, i]
        avg = np.average(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg
    return dict



def plot_v_lightcurve_with_fit(SN_dict, sampler):
    param_dict = get_param_results_dict(sampler)
    m0 = param_dict['m0']
    a0 = param_dict['a0']
    tPT = param_dict['tPT']
    w0 = param_dict['w0']

    # print('m0:', m0, 'a0:', a0, 'tPT:', tPT, 'w0:', w0)
    p0 = 0.013174860582548461
    data = get_LCO_V_df(SN_dict)
    data = data.loc[data['t_from_discovery'] < 190]
    x = data['t_from_discovery']
    y = data['mag']
    dy = data['dmag']
    y_fit = m0 + p0 * x - a0 / (1 + np.exp((x - tPT) / w0))
    plt.figure()
    plt.errorbar(x, y, yerr=dy, marker='o', linestyle='None')
    plt.plot(x, y_fit)
    plt.gca().invert_yaxis()
