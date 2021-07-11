import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec
import mcmc_snec_noCSM

'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

SN_name = 'SN2017eaw'
extra_tail_days = False  # has to be False or an integer

distances = pd.read_csv(os.path.join('results', 'distances.csv'))
distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
sigma_S = 2 * distance_err / distance

Mzams_range = [9.0, 11.0, 13.0, 15.0, 17]
Ni_range = [0.02, 0.07, 0.12]
E_final_range = [0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]
R_range = [0, 500, 1000, 2000]
K_range = [0, 10, 30, 60]
S_range = [1.0-sigma_S, 1.0+sigma_S]
T_range = [-10, 2]

# nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

n_walkers = 50
n_steps = 800
n_params = 6
burn_in = 1200

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_twostep_'+SN_name)
Path(res_dir).mkdir(parents=True, exist_ok=True)

run_param_df = pd.DataFrame.from_dict({'SN_name': SN_name,
                             'Mzams_range': str(Mzams_range), 'Ni_range': str(Ni_range),
                             'E_final_range': str(E_final_range), 'R_range': str(R_range), 'K_range': str(K_range),
                             'Mix_range': str(Mix_range), 'Scaling_range': str(S_range), 'T_range': str(T_range),
                             'n_walkers': n_walkers, 'n_steps_per_stage': n_steps, 'n_params': n_params,
                             'burn_in': burn_in,
                             'time': time_now}, orient='index')
run_param_df.to_csv(os.path.join(res_dir, 'run_parameters.csv'))

m_Solar = 1.989 * (10 ** 33)  # gram

parameter_ranges = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'Mix' : Mix_range, 'S' : S_range, 'T' : T_range}

# import SN bolometric lum data
data_lum= pd.read_csv(os.path.join('results', SN_name+'_martinez.csv'))
# take only days between 30 and 200 (like the Martinez paper)
if extra_tail_days is not False:
    data_lum = data_lum.loc[data_lum['t_from_discovery'] <= float(200 + extra_tail_days - 10)]
else:
    data_lum = data_lum.loc[data_lum['t_from_discovery'] <= 200]
data_lum_no_early = data_lum.loc[data_lum['t_from_discovery'] > 30]
# data_lum= pd.read_csv(os.path.join('results', 'blackbody_results_'+SN_name+'_BVgri.csv'))
# convert watt to erg/s
# data_lum['Lum'] = data_lum['Lum'] * 10**7
# data_lum['dLum0'] = data_lum['dLum0'] * 10**7
# data_lum['dLum1'] = data_lum['dLum1'] * 10**7


# import SN photospheric velocities lum data
data_veloc = pd.read_csv(os.path.join('results', SN_name+'_expansion_velocities.csv'),
                      usecols=['t_from_discovery', 'line', 'absorption_mean_velocity','absorption_std_velocity'])
data_veloc.rename({'absorption_mean_velocity':'veloc', 'absorption_std_velocity':'dveloc'}, axis='columns', inplace=True)
data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
# remove first point which seems like an artifact
# data_veloc = data_veloc.loc[data_veloc['t_from_discovery'] > 20]


'''
# running code #
'''

SN_data_all = {'veloc': data_veloc, 'lum': data_lum}
SN_data_all_no_early = {'veloc': data_veloc, 'lum': data_lum_no_early}


"""
first step: fit to Mzams, Ni, E, Mix, S, T, without CSM. With gaussian prior for S, uniform for all others.
"""


sampler_first_step = mcmc_snec_noCSM.emcee_fit_params(SN_data_all_no_early, n_walkers, n_steps, parameter_ranges,
                                     'lum_veloc_Tthresh_normalized', extend_tail=extra_tail_days)


# adding zero CSM columns to the sampler chain of the first step

sampler_first_step_chain_flat = sampler_first_step.get_chain(flat=True)
sampler_first_step_chain_flat_CSM = np.zeros((n_steps*n_walkers, 8))
sampler_first_step_chain_flat_CSM[:,0] = sampler_first_step_chain_flat[:,0]
sampler_first_step_chain_flat_CSM[:,1] = sampler_first_step_chain_flat[:,1]
sampler_first_step_chain_flat_CSM[:,2] = sampler_first_step_chain_flat[:,2]
sampler_first_step_chain_flat_CSM[:,5] = sampler_first_step_chain_flat[:,3]
sampler_first_step_chain_flat_CSM[:,6] = sampler_first_step_chain_flat[:,4]
sampler_first_step_chain_flat_CSM[:,7] = sampler_first_step_chain_flat[:,5]

sampler_first_step_chain = sampler_first_step.chain
sampler_first_step_chain_CSM = np.zeros((n_walkers, n_steps, 8))
sampler_first_step_chain_CSM[:,:,0] = sampler_first_step_chain[:,:,0]
sampler_first_step_chain_CSM[:,:,1] = sampler_first_step_chain[:,:,1]
sampler_first_step_chain_CSM[:,:,2] = sampler_first_step_chain[:,:,2]
sampler_first_step_chain_CSM[:,:,5] = sampler_first_step_chain[:,:,3]
sampler_first_step_chain_CSM[:,:,6] = sampler_first_step_chain[:,:,4]
sampler_first_step_chain_CSM[:,:,7] = sampler_first_step_chain[:,:,5]

"""
second step: fit to K and R, with priors according to results from step one.
"""

flat_sampler_no_burnin = pd.DataFrame(sampler_first_step_chain_flat_CSM[n_walkers*burn_in:,:])
flat_sampler_no_burnin.columns = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']


nonuniform_priors = {}
for param in ['Mzams', 'Ni', 'E', 'Mix', 'S', 'T']:
    nonuniform_priors[param] = {'polynomial': mcmc_snec.polyfit_to_distribution(flat_sampler_no_burnin[param], res_dir)}

parameter_ranges = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range, 'S' : S_range, 'T' : T_range}

last_walkers_step_one = sampler_first_step_chain_flat_CSM[-n_walkers:, :]
# change R and K to initial guesses of uniform dist within range
last_walkers_step_one[:, 3] = np.random.rand(n_walkers) * \
                              (parameter_ranges['R'][-1] - parameter_ranges['R'][0]) + parameter_ranges['R'][0]
last_walkers_step_one[:, 4] = np.random.rand(n_walkers) * \
                              (parameter_ranges['K'][-1] - parameter_ranges['K'][0]) + parameter_ranges['K'][0]

sampler_second_step = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                                 'lum_veloc_Tthresh_normalized',
                                                 nonuniform_priors, init_guesses=last_walkers_step_one,
                                                 extend_tail=extra_tail_days)


sampler_chain = np.concatenate((sampler_first_step_chain_CSM, sampler_second_step.chain), axis=1)
param_dict = mcmc_snec.get_param_results_dict(sampler_chain, parameter_ranges, n_steps*2-1, res_dir)

sampler_chain_flat = np.concatenate((
    sampler_first_step_chain_flat_CSM,
    sampler_second_step.get_chain(flat=True)), axis=0)

pd.DataFrame(sampler_chain_flat).to_csv(os.path.join(res_dir, 'flat_sampler.csv'))


mcmc_snec.plot_lightcurve_with_fit(sampler_chain, SN_data_all, parameter_ranges,
                                        'lum_veloc_Tthresh_normalized', res_dir, n_walkers, SN_name, n_steps*2-1, extend_tail=extra_tail_days)
mcmc_snec.plot_lightcurve_with_fit(sampler_chain, SN_data_all, parameter_ranges,
                                        'lum_veloc_Tthresh_normalized', res_dir, n_walkers, SN_name, 0, extend_tail=extra_tail_days)

mcmc_snec.chain_plots(sampler_chain, parameter_ranges, res_dir, burn_in, first_stage_steps=n_steps)

mcmc_snec.corner_plot(sampler_chain_flat[n_walkers*burn_in:, :], parameter_ranges, res_dir)


