import emcee
import numpy as np
from matplotlib import pyplot as plt

def log_prior(theta):
    # theta is a vector containing a specific set of parameters theta = [m0, a0, tPT, w0]
    m_Ni_min = 0
    m_Ni_max = 10
    if m_Ni_min < theta[0] < m_Ni_max:
        prob_m_Ni = 1. / theta[0]
    else:
        prob_m_Ni = 0.
    prob_total = np.log(prob_m_Ni)
    return prob_total

m_Solar = 1.989 * (10 ** 33)  # gram
print('solar', m_Solar)

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
    n_steps = 300
    n_params = 1
    args = [bolometric_time_vec, bolometric_mag_vec, bolometric_dmag_vec]

    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior, args=args)

    m_Ni_min = 0
    m_Ni_max = 10

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


def plot_v_lightcurve_with_fit(blackbody_results, sampler):
    Ni = np.average(sampler.chain[:, -1, 0])

    Ni_std = np.std(sampler.chain[:, -1, 0])

    print(Ni)
    print(Ni_std)

    data = blackbody_results
    x = data['t_from_discovery']
    y = data['Lum'] / 10 ** 43
    dy = data['dLum0'] / 10 ** 43
    y_fit = Ni * (6.45 * np.exp(-x / 8.8) + 1.45 * np.exp(-x / 111.3))
    plt.figure()
    plt.errorbar(x, y, yerr=dy, marker='o', linestyle='None', label='bolometric luminosity')
    plt.ylabel('bolometric luminosity (10^43 erg/s)')
    plt.xlabel('time after discovery (days)')
    plt.plot(x, y_fit, label='Ni mass fit')
    plt.legend()
    return y_fit

