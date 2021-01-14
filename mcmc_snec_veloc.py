import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec

'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

SN_name = 'SN2012ec'

Mzams_range = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
Ni_range = [0.02, 0.07, 0.12, 0.17]
E_final_range = [0.7, 0.9, 1.1, 1.3, 1.5]
Mix_range = [5.0, 8.0]
R_range = [0.1, 0.1]
K_range = [0.1, 0.1]
S_range = [0.8, 1.2]
T_range = [-15, 15]

n_walkers = 30
n_steps = 300
n_params = 8
burn_in = 100

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_veloc')
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


# import SN photospheric velocities data
data_veloc = pd.read_csv(os.path.join('results', SN_name+'_expansion_velocities.csv'),
                      usecols=['t_from_discovery', 'line', 'absorption_mean_velocity','absorption_std_velocity'])
data_veloc.rename({'absorption_mean_velocity':'veloc', 'absorption_std_velocity':'dveloc'}, axis='columns', inplace=True)
data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
# remove first point which seems like an artifact
data_veloc = data_veloc.loc[data_veloc['t_from_discovery'] > 20]
print(data_veloc)



'''
# running code #
'''

SN_data_all = {'veloc': data_veloc}
SN_data_all_w_early = {'veloc': data_veloc}

sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges, 'veloc')

flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")

mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all_w_early, parameter_ranges,
                                        'veloc', res_dir, n_walkers, SN_name, n_steps-1)
mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all_w_early, parameter_ranges,
                                        'veloc', res_dir, n_walkers, SN_name, 0)


mcmc_snec.chain_plots(sampler, parameter_ranges, res_dir)

mcmc_snec.corner_plot(sampler, burn_in, parameter_ranges, res_dir)

print(sampler.chain.shape)

