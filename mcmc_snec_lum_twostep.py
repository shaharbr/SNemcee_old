import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec
from matplotlib import pyplot as plt
from numpy import trapz



time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_lum_twostep')
Path(res_dir).mkdir(parents=True, exist_ok=True)

def polyfit_to_distribution(array_walker_results):
    name = array_walker_results.name
    counts, bin_edges = np.histogram(array_walker_results, bins=20)
    fig, ax = plt.subplots()
    ax.hist(array_walker_results, density=True, bins=20)
    bin_widths = np.diff(bin_edges)
    x = bin_edges[:-1] + (bin_widths / 2)
    y = counts
    area = trapz(y, dx=bin_widths[0])
    y = y / area
    polymod = np.poly1d(np.polyfit(x, y, deg=12))
    ax.plot(x, y, color='black')
    dense_x = np.arange(np.min(x), np.max(x), (np.max(x)-np.min(x))/50)
    ax.plot(dense_x, [polymod(i) for i in dense_x], color='orange')
    fig.savefig(os.path.join(res_dir, str(name)+'_first_step_dist.png'))
    return polymod



'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

SN_name = 'SN2012ec'

distances = pd.read_csv(os.path.join('results', 'distances.csv'))
distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
sigma_S = 2 * distance_err / (distance)


Mzams_range = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
Ni_range = [0.02, 0.07, 0.12]
E_final_range = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]
Mix_range = [2.0, 5.0, 8.0]
R_range = [0, 0]
K_range = [0, 0]
S_range = [1.0-sigma_S, 1.0+sigma_S]
T_range = [-10, 2]

nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

n_walkers = 70
n_steps = 400
n_params = 8
burn_in = 150


run_param_df = pd.DataFrame.from_dict({'SN_name': SN_name,
                             'Mzams_range': str(Mzams_range), 'Ni_range': str(Ni_range),
                             'E_final_range': str(E_final_range), 'R_range': str(R_range), 'K_range': str(K_range),
                             'Mix_range': str(Mix_range), 'Scaling_range': str(S_range), 'T_range': str(T_range),
                             'n_walkers': n_walkers, 'n_steps': n_steps, 'n_params': n_params,
                             'burn_in': burn_in,
                             'time': time_now}, orient='index')
run_param_df.to_csv(os.path.join(res_dir, 'run_parameters.csv'))

m_Solar = 1.989 * (10 ** 33)  # gram

parameter_ranges = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range, 'S' : S_range, 'T' : T_range}

# import SN bolometric lum data
data_lum_w_early= pd.read_csv(os.path.join('results', SN_name+'_martinez.csv'))
# convert watt to erg/s
# data_lum_w_early['Lum'] = data_lum_w_early['Lum'] * 10**7
# data_lum_w_early['dLum0'] = data_lum_w_early['dLum0'] * 10**7
# data_lum_w_early['dLum1'] = data_lum_w_early['dLum1'] * 10**7
# take only days between 30 and 200 (like the Martinez paper)
data_lum_w_early = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] < 250]

# replicate the last point x times to artificially increase weight of fitting to last point (10-folds)
times_to_amplify = 1
if times_to_amplify > 1:
    last_row = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] > 350]
    last_row_repeats = pd.concat([last_row]*(times_to_amplify-1), ignore_index=True)
    data_lum_w_early= pd.concat([data_lum_w_early, last_row_repeats], ignore_index=True)

data_lum = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] > 30]


'''
# running code #
'''

SN_data_all = {'lum': data_lum}
SN_data_all_w_early = {'lum': data_lum_w_early}



"""
first step: fit to Mzams, Ni, E, Mix, S, T, without CSM. With gaussian prior for S, uniform for all others.
"""


sampler_first_step = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges, 'lum', nonuniform_priors)


"""
second step: fit to K and R, with priors according to results from step one.
"""

flat_sampler_no_burnin = pd.DataFrame(sampler_first_step.get_chain(discard=burn_in, flat=True))
flat_sampler_no_burnin.columns = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
for param in ['Mzams', 'Ni', 'E', 'Mix', 'S', 'T']:
    nonuniform_priors[param] = {'polynomial': polyfit_to_distribution(flat_sampler_no_burnin[param])}

parameter_ranges['Mzams'] = [9.0, 11.0, 13.0, 15.0, 17.0]
parameter_ranges['Ni'] = [0.02, 0.12]
parameter_ranges['E'] = [0.5, 0.9, 1.3, 1.7]
parameter_ranges['Mix'] = [2.0, 8.0]
parameter_ranges['R'] = [1000, 2000]
parameter_ranges['K'] = [0, 20, 50, 100]

last_walkers_step_one = sampler_first_step.get_chain(flat=True)[-(n_walkers+1):-1, :]

# change R and K to initial guesses of uniform dist within range
last_walkers_step_one[:, 3] = np.random.rand(n_walkers) * \
                              (parameter_ranges['R'][-1] - parameter_ranges['R'][0]) + parameter_ranges['R'][0]
last_walkers_step_one[:, 4] = np.random.rand(n_walkers) * \
                              (parameter_ranges['K'][-1] - parameter_ranges['K'][0]) + parameter_ranges['K'][0]



sampler_second_step = mcmc_snec.emcee_fit_params(SN_data_all_w_early, n_walkers, n_steps, parameter_ranges, 'lum',
                                                 nonuniform_priors, init_guesses=last_walkers_step_one)
sampler_chain = np.concatenate((sampler_first_step.chain, sampler_second_step.chain), axis=1)
param_dict = mcmc_snec.get_param_results_dict(sampler_chain, parameter_ranges, n_steps*2-1, res_dir)

sampler_chain_flat = np.concatenate((
    sampler_first_step.get_chain(flat=True),
    sampler_second_step.get_chain(flat=True)), axis=0)


pd.DataFrame(sampler_chain_flat).to_csv(os.path.join(res_dir, 'flat_sampler.csv'))

mcmc_snec.plot_lightcurve_with_fit(sampler_chain, SN_data_all_w_early, parameter_ranges,
                                        'lum', res_dir, n_walkers, SN_name, n_steps*2-1)
mcmc_snec.plot_lightcurve_with_fit(sampler_chain, SN_data_all_w_early, parameter_ranges,
                                        'lum', res_dir, n_walkers, SN_name, 0)

mcmc_snec.chain_plots(sampler_chain, parameter_ranges, res_dir, burn_in, first_stage_steps=n_steps)

mcmc_snec.corner_plot(sampler_chain_flat[n_walkers*burn_in-1:-1, :], parameter_ranges, res_dir)

