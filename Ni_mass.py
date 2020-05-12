import emcee
import numpy as np
from matplotlib import pyplot as plt


# TODO in scipy too, least suare fitting
# TODO check about the two outlier bolometric points, whats going on in the filters


plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



def log_prior(theta):
    # theta is a vector containing a specific set of parameters theta = [m0, a0, tPT, w0]
    m_Ni_min = 0.0
    m_Ni_max = 0.3
    if m_Ni_min < theta[0] < m_Ni_max:
        prob_m_Ni = 1. / theta[0]
    else:
        prob_m_Ni = 0.
    prob_total = np.log(prob_m_Ni)
    return prob_total

m_Solar = 1.989 * (10 ** 33)  # gram

def log_likelihood(theta, data_x, data_y, data_dy):
    [m_Ni] = theta
    y_fit = m_Ni * (6.45 * np.exp(-data_x/8.8) + 1.45 * np.exp(-data_x/111.3))* 10 ** 43
    chi2 = (data_y - y_fit) ** 2. / (2. * data_dy ** 2.) + np.log(data_y)
    ln_like = -np.sum(chi2)
    return ln_like


def log_posterior(theta, data_x, data_y, data_dy):
    ln_post = log_prior(theta) + log_likelihood(theta, data_x, data_y, data_dy)
    return ln_post



def emcee_fit_params(bolometric_time_vec, bolometric_mag_vec, bolometric_dmag_vec):
    n_walkers = 50
    n_steps = 100
    n_params = 1
    args = [bolometric_time_vec, bolometric_mag_vec, bolometric_dmag_vec]

    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=args)

    m_Ni_min = 0.0
    m_Ni_max = 0.3

    m_Ni_random = np.random.rand(n_walkers) * (m_Ni_max - m_Ni_min) + m_Ni_min
    initial_guesses = np.array([m_Ni_random]).T

    sampler.run_mcmc(initial_guesses, n_steps)

    return sampler


def SN_lightcurve_params(blackbody_results):
    # take only Ni tail
    blackbody_results = blackbody_results.loc[blackbody_results['t_from_discovery'] > 140]
    bolometric_time_vec = blackbody_results['t_from_discovery']
    bolometric_mag_vec = blackbody_results['Lum']
    bolometric_dmag_vec = blackbody_results['dLum0']
    sampler = emcee_fit_params(bolometric_time_vec, bolometric_mag_vec, bolometric_dmag_vec)
    return sampler


def chain_plots(chain, **kwargs):
    plt.figure()
    plt.plot(chain[:, :, 0].T, **kwargs)
    plt.xlabel('Step Number')
    plt.ylabel('Ni')

def get_Ni_results_dict(sampler):
    last_results = sampler.chain[:, -50:, 0]
    avg = np.average(last_results)
    sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
    sigma_lower = avg - sigma_lower
    sigma_upper = sigma_upper - avg
    dict = {'Ni_mcmc': avg, 'Ni_mcmc_16perc': sigma_lower, 'Ni_mcmc_84perc': sigma_upper}
    return dict




def plot_v_lightcurve_with_fit(blackbody_results, sampler):
    Ni = np.average(sampler.chain[:, -1, 0])

    Ni_std = np.std(sampler.chain[:, -1, 0])

    # print(Ni)
    # print(Ni_std)
    data = blackbody_results
    x = data['t_from_discovery']
    y = data['Lum'] / 10 ** 43
    dy = data['dLum0'] / 10 ** 43
    y_fit = Ni * (6.45 * np.exp(-x / 8.8) + 1.45 * np.exp(-x / 111.3))
    fig, ax = plt.subplots(1, figsize=(5, 5))
    ax.errorbar(x, y, yerr=dy, marker='o', linestyle='None', label='bolometric luminosity')
    ax.set_ylabel('bolometric luminosity (10^43 erg/s)')
    ax.set_xlabel('Rest-frame days from discovery')
    ax.plot(x, y_fit, label='Ni fit')
    ax.set_yscale('log')
    ax.legend()
    fig.savefig(r'figures\SN2018hmx_Ni_mass_luminosity_fit_to_bolometric' + '.png')
    fig.savefig(r'figures\SN2018hmx_Ni_mass_luminosity_fit_to_bolometric' + '.svg')
    return y_fit



# def log_slope_ratio(slope1, slope2):
#     return np.power(10, slope1 - slope2)

# def calc_87A_slope_ratio(day, SN_fit, SN87A_fit):
#     val = 0.075 * SN_fitnp.power(10, SN_line[1] + day * SN_line[0]) / np.power(10, SN87A_line[1] + day * SN87A_line[0]))
#     return val


def fit_to_log_slope(SN_bolometric_df):
    tail = SN_bolometric_df.loc[SN_bolometric_df['t_from_discovery'] > 132]
    x = tail['t_from_discovery']
    y = np.log10(tail['Lum'])
    if 'dLum0' in SN_bolometric_df.columns:
        yerr = np.log10(tail['dLum0'])
        weights = list(1 / yerr)
    else:
        weights = None
    line, cov = np.polyfit(x, y, deg=1, w=weights, cov=True)

    x_fit = np.arange(132, 400, 10)
    y_fit_log = line[1] + x_fit * line[0]

    sigma_log = np.sqrt(np.diag(cov))
    # print('line', line)
    # print('cov', cov, sigma_log)
    sigma_log = np.sqrt((sigma_log[0] * x_fit) ** 2 + sigma_log[1] ** 2)

    y_fit = np.power(10, y_fit_log)
    y_fit_sigma = sigma_log * np.log(10) * np.power(10, y_fit_log)
    # print('yfit')
    # print(y_fit_log)
    # print(y_fit)

    # print('sigma')
    # print(sigma_log)
    # print(y_fit_sigma)
    return x_fit, y_fit, y_fit_sigma



def calc_87A_slope_ratio(SN_fit, SN87A_fit):
    Ni = 0.075 * SN_fit / SN87A_fit
    return Ni

def calc_87A_slope_ratio_sigma(SN_fit, SN87A_fit, SN_fit_sigma, SN87A_fit_sigma):
    Ni_sigma = np.sqrt((SN_fit_sigma * 0.075 / SN87A_fit) ** 2 + (SN87A_fit_sigma * 0.075 * SN_fit / (SN87A_fit ** 2)) ** 2)
    return Ni_sigma


def Ni_by_87A_slope(SN_fit, SN_fit_sigma, SN87A_fit, SN87A_fit_sigma):
    Ni_mass = np.average([calc_87A_slope_ratio(SN_fit[i], SN87A_fit[i])
                          for i in range(1)])
    Ni_sigma = np.average([calc_87A_slope_ratio_sigma(SN_fit[i], SN87A_fit[i], SN_fit_sigma[i], SN87A_fit_sigma[i])
                           for i in range(1)])
    return Ni_mass, Ni_sigma


# def calc_87A_slope(SN1987A_bolometric):
    # fit
    # SN87A_tail = SN1987A_bolometric.loc[SN1987A_bolometric['t_from_discovery'] > 132]
    # x_87A = SN87A_tail['t_from_discovery']
    # y_87A = np.log10(SN87A_tail['Lum'])
    # SN87A_line = np.polyfit(x_87A, y_87A, deg=1)
    # SN87A_line, SN87A_cov = np.polyfit(x_87A, y_87A, deg=1, cov=True)
    # SN87A_sigma = np.sqrt(np.diag(SN87A_cov))
    # SN87A_sigma = SN87A_sigma[0]
    # print('sig', SN87A_sigma)
    # return SN87A_line, SN87A_sigma


