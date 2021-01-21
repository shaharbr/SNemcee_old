import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec

'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

SN_name = 'SN2017eaw'

distances = pd.read_csv(os.path.join('results', 'distances.csv'))
distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
sigma_S = 2 * distance_err / (distance)


Mzams_range = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
Ni_range = [0.02, 0.07, 0.12, 0.17]
E_final_range = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]
Mix_range = [2.0, 5.0, 8.0]
R_range = [0, 0]
K_range = [0, 0]

# Mzams_range = [9.0, 11.0, 13.0, 15.0, 17.0]
# Ni_range = [0.02, 0.12]
# E_final_range = [0.5, 0.9, 1.3, 1.7]
# Mix_range = [2.0, 8.0]
# R_range = [1000, 2000]
# K_range = [0, 50, 100]
S_range = [1.0-sigma_S, 1.0+sigma_S]
T_range = [-10, 2]

nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

n_walkers = 70
n_steps = 800
n_params = 8
burn_in = 150

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_lum_'+SN_name)
Path(res_dir).mkdir(parents=True, exist_ok=True)

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
data_lum_w_early = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] < 200]

# replicate the last point x times to artificially increase weight of fitting to last point (10-folds)
times_to_amplify = 1
if times_to_amplify > 1:
    last_row = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] > 350]
    last_row_repeats = pd.concat([last_row]*(times_to_amplify-1), ignore_index=True)
    data_lum_w_early= pd.concat([data_lum_w_early, last_row_repeats], ignore_index=True)

# data_lum = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] > 30]
data_lum = data_lum_w_early

'''
# running code #
'''

SN_data_all = {'lum': data_lum}
SN_data_all_w_early = {'lum': data_lum_w_early}


sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges, 'lum', nonuniform_priors)

flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")

mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all_w_early, parameter_ranges,
                                        'lum', res_dir, n_walkers, SN_name, n_steps-1)
mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all_w_early, parameter_ranges,
                                        'lum', res_dir, n_walkers, SN_name, 0)

mcmc_snec.chain_plots(sampler.chain, parameter_ranges, res_dir, burn_in)

mcmc_snec.corner_plot(sampler.get_chain(flat=True)[n_walkers*burn_in:-1, :], parameter_ranges, res_dir)

print(sampler.chain.shape)

