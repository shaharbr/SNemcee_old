import emcee
import numpy as np
from matplotlib import pyplot as plt


def log_prior(theta):
    # theta is a vector containing a specific set of parameters theta = [m0, a0, tPT, w0]
    tPT_min = 50
    tPT_max = 150
    a0_min = 0
    a0_max = 6
    w0_min = 0
    w0_max = 100
    m0_min = 10
    m0_max = 25
    if m0_min < theta[0] < m0_max:
        prob_m0 = 1. / theta[0]
    else:
        prob_m0 = 0.
    if a0_min < theta[1] < a0_max:
        prob_a0 = 1. / theta[1]
    else:
        prob_a0 = 0.
    if tPT_min < theta[2] < tPT_max:
        prob_tPT = 1. / theta[2]
    else:
        prob_tPT = 0.
    if w0_min < theta[3] < w0_max:
        prob_w0 = 1. / theta[3]
    else:
        prob_w0 = 0.
    prob_total = np.log(prob_m0 * prob_a0 * prob_tPT * prob_tPT * prob_w0)
    return prob_total

def log_likelihood(theta, data_x, data_y, data_dy):
    [m0, a0, tPT, w0] = theta
    p0 = 0.013174860582548461
    y_fit = m0 + p0 * data_x - a0 / (1 + np.exp((data_x - tPT) / w0))
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



def calc_Vmag_slope_param(data, time_range):
    first_day = time_range[0]
    last_day = time_range[1]
    data = data.loc[(data['t_from_discovery'] >= first_day) &
                    (data['t_from_discovery'] <= last_day)]
    x = list(data['t_from_discovery'])
    y = list(data['mag'])
    weights = list(1/data['dmag'])
    regression_params, cov = np.polyfit(x, y, deg=1, w=weights, cov=True)
    # print('polyfit', regression_params, cov)
    sigma = np.sqrt(np.diag(cov))
    # print('sigma', sigma)
    return regression_params, sigma

def calc_s50V(SN_dict, get_all_params=False):
    data = get_LCO_V_df(SN_dict)
    peak_mag = np.min(data.loc[data['t_from_discovery'] < 100, 'mag'])
    first_day = float(data.loc[data['mag'] == peak_mag, 't_from_discovery'])
    last_day = first_day + 97
    V50_regression_params, sigma = calc_Vmag_slope_param(data, time_range=[first_day, last_day])
    if get_all_params:
        s50V = V50_regression_params
    else:
        s50V = V50_regression_params[0] * 50
    sigma = sigma[0]
    print('s50V:', s50V)
    return s50V, sigma

def calc_p0(SN_dict, time_range=False, get_all_params=False):
    data = get_LCO_V_df(SN_dict)
    if time_range:
        first_day = time_range[0]
        last_day = time_range[1]
    else:
        last_day = np.max(data['t_from_discovery'])
        first_day = last_day - 30
    p0_regression_params, sigma = calc_Vmag_slope_param(data, time_range=[first_day, last_day])
    if get_all_params:
        p0 = p0_regression_params
    else:
        p0 = p0_regression_params[0]
    sigma = sigma[0]
    print('p0:', p0)
    return p0, sigma


def plot_v_lightcurve_with_slope(SN_dict, param):
    data = get_LCO_V_df(SN_dict)
    data_x = data['t_from_discovery']
    data_y = data['mag']
    data_dy = data['dmag']
    if param == 'p0':
        time_range = [140, 190]
        regression_params = calc_p0(SN_dict, time_range, get_all_params=True)
    elif param == 's50V':
        time_range = [15, 65]
        regression_params = calc_s50V(SN_dict, get_all_params=True)
    slope_x = time_range
    slope_y = [regression_params[0] * t + regression_params[1] for t in slope_x]

    plt.figure()
    plt.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None')
    plt.plot(slope_x, slope_y, color='purple')
    plt.gca().invert_yaxis()
    # TODO fix


def calc_Vmax(SN_dict):
    V_df = SN_dict['lightcurve']['Las Cumbres']['df']
    Vmag = V_df.loc[V_df['filter'] == 'V', 'abs_mag']
    Vmax = np.min(Vmag)
    Vmax_err = np.min(V_df.loc[(V_df['filter'] == 'V') & (V_df['abs_mag'] == Vmax), 'dmag'])
    print('Vmax: ', Vmax)
    print('Vmax_err: ', Vmax_err)
    return Vmax, Vmax_err



