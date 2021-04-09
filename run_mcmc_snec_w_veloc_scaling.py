import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec_w_veloc_scaling as mcmc_snec

'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

SN_name = 'SN2012aw'
extra_tail_days = 60  # has to be False or an integer

distances = pd.read_csv(os.path.join('results', 'distances.csv'))
distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
sigma_S = 2 * distance_err / distance

Mzams_range = [9.0, 11.0, 13.0, 15.0]
Ni_range = [0.02, 0.07, 0.12]
E_final_range = [0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]
R_range = [500, 1000, 2000]
K_range = [10, 30, 60]
S_range = [1.0-sigma_S, 1.0+sigma_S]
T_range = [-10, 2]
Sv_range = [0.5, 2.0]

# nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

n_walkers = 100
n_steps = 1500
n_params = 9
burn_in = 500

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_lum_veloc_w_Sv_'+SN_name)
Path(res_dir).mkdir(parents=True, exist_ok=True)

run_param_df = pd.DataFrame.from_dict({'SN_name': SN_name,
                             'Mzams_range': str(Mzams_range), 'Ni_range': str(Ni_range),
                             'E_final_range': str(E_final_range), 'R_range': str(R_range), 'K_range': str(K_range),
                             'Mix_range': str(Mix_range), 'Scaling_range': str(S_range), 'T_range': str(T_range),
                             'velocity_scaling_range': str(Sv_range),
                             'n_walkers': n_walkers, 'n_steps': n_steps, 'n_params': n_params,
                             'burn_in': burn_in,
                             'time': time_now}, orient='index')
run_param_df.to_csv(os.path.join(res_dir, 'run_parameters.csv'))

m_Solar = 1.989 * (10 ** 33)  # gram


parameter_ranges = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range, 'S' : S_range, 'T' : T_range,
                    'Sv' : Sv_range}

# import SN bolometric lum data
data_lum= pd.read_csv(os.path.join('results', SN_name+'_martinez.csv'))
# take only days between 30 and 200 (like the Martinez paper)
if extra_tail_days is not False:
    data_lum = data_lum.loc[data_lum['t_from_discovery'] <= float(200 + extra_tail_days - 10)]
else:
    data_lum = data_lum.loc[data_lum['t_from_discovery'] <= 200]
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

sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                     'lum_veloc_Tthresh_normalized', extend_tail=extra_tail_days)

flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")

mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all, parameter_ranges,
                                        'lum_veloc_Tthresh_normalized', res_dir, n_walkers, SN_name, n_steps-1, extend_tail=extra_tail_days)
mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all, parameter_ranges,
                                        'lum_veloc_Tthresh_normalized', res_dir, n_walkers, SN_name, 0, extend_tail=extra_tail_days)

mcmc_snec.chain_plots(sampler.chain, parameter_ranges, res_dir, burn_in)

mcmc_snec.corner_plot(sampler.get_chain(flat=True)[n_walkers*burn_in:-1, :], parameter_ranges, res_dir)

print(sampler.chain.shape)
