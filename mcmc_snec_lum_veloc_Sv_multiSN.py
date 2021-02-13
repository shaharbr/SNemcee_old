import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec_Sv_Tthresh_multiSN as mcmc_snec




'''
parts of this code are based on code by Griffin Hosseinzadeh
'''

# SN_names = ['SN2004a', 'SN2005cs', 'SN2012aw', 'SN2012ec', 'SN2017eaw', 'SN2018aoq', 'SN2004et_d5.9']
# SN_names = ['SN2004a', 'SN2012aw', 'SN2012ec', 'SN2017eaw', 'SN2004et_d5.9']
SN_names = ['SN2004a', 'SN2012aw']
low_luminosity_SNe = ['SN2005cs', 'SN2018aoq']
num_SNe = len(SN_names)

Mzams_range = [9.0, 11.0, 13.0, 15.0, 17.0]
Ni_range = [0.02, 0.12]
E_final_range = [0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]

# Mzams_range = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
# Ni_range = [0.02, 0.07, 0.12]
# E_final_range = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]
# Mix_range = [2.0, 5.0, 8.0]

Mzams_range_low = [9.0, 10.0, 11.0, 0.12]
Ni_range_low = [0.001, 0.02, 0.07, 0.12, 0.17]
E_final_range_low = [0.1, 0.3, 0.5, 0.7, 0.9]
Mix_range_low = [2.0, 5.0, 8.0]

R_range = [1000, 2000]
K_range = [0, 20, 50, 100]
S_range = [0.5, 2.0]
T_range = [-10, 2]
Sv_range = [0.1, 2.0]

param_high = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range, 'S' : S_range, 'T' : T_range,
                    'Sv' : Sv_range}

param_low = {'Mzams' : Mzams_range_low, 'Ni' : Ni_range_low, 'E' : E_final_range_low,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range_low, 'S' : S_range, 'T' : T_range,
                    'Sv' : Sv_range}

n_steps = 100
n_params = 8 * num_SNe + 1
n_walkers = (n_params + 1) * 2
burn_in = 50

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_lum_veloc_findingSv_allSNe')
Path(res_dir).mkdir(parents=True, exist_ok=True)

data_lum_list = []
data_lum_w_early_list = []
data_veloc_list = []
parameter_ranges = []

for SN_name in SN_names:
    # import SN bolometric lum data
    data_lum_w_early = pd.read_csv(os.path.join('results', SN_name + '_martinez.csv'))
    # take only days between 30 and 200 (like the Martinez paper)
    data_lum_w_early = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] < 200]
    data_lum = data_lum_w_early.loc[data_lum_w_early['t_from_discovery'] > 30]
    # import SN photospheric velocities lum data
    data_veloc = pd.read_csv(os.path.join('results', SN_name + '_expansion_velocities.csv'),
                             usecols=['t_from_discovery', 'line', 'absorption_mean_velocity',
                                      'absorption_std_velocity'])
    data_veloc.rename({'absorption_mean_velocity': 'veloc', 'absorption_std_velocity': 'dveloc'}, axis='columns',
                      inplace=True)
    data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
    data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
    # remove first point which seems like an artifact
    data_veloc = data_veloc.loc[data_veloc['t_from_discovery'] > 20]
    # add to list of all SNe data
    data_lum_list.append(data_lum)
    data_lum_w_early_list.append(data_lum_w_early)
    data_veloc_list.append(data_veloc)
    sigma_S = mcmc_snec.sigma_S(SN_name)
    if SN_name in low_luminosity_SNe:
        param_low['S'] = [1.0 - sigma_S, 1.0 + sigma_S]
        parameter_ranges.append(param_low)
    else:
        param_high['S'] = [1.0 - sigma_S, 1.0 + sigma_S]
        parameter_ranges.append(param_high)

# print(data_lum_list)
# print(data_lum_w_early_list)
# print(data_veloc_list)
# print(parameter_ranges)


run_param_df = pd.DataFrame.from_dict({'SN_name': str(SN_names),
                             'Mzams_range_highlum': str(Mzams_range),
                             'Ni_range_highlum': str(Ni_range),
                             'E_final_range_highlum': str(E_final_range),
                             'Mix_range_highlum': str(Mix_range),
                             'Mzams_range_lowlum': str(Mzams_range),
                             'Ni_range_lowlum': str(Ni_range_low),
                             'E_final_range_lowlum': str(E_final_range_low),
                             'Mix_range_lowlum': str(Mix_range_low),
                             'R_range': str(R_range), 'K_range': str(K_range),
                             'Scaling_range': str(S_range), 'T_range': str(T_range),
                             'velocity_scaling_range': str(Sv_range),
                             'n_walkers': n_walkers, 'n_steps': n_steps, 'n_params': n_params,
                             'burn_in': burn_in,
                             'time': time_now}, orient='index')
run_param_df.to_csv(os.path.join(res_dir, 'run_parameters.csv'))

'''
# running code #
'''


SN_data_all = [{'veloc': data_veloc_list[i], 'lum': data_lum_list[i]} for i in range(num_SNe)]
SN_data_all_w_early = [{'veloc': data_veloc_list[i], 'lum': data_lum_w_early_list[i]} for i in range(num_SNe)]

sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges, SN_names, 'lum_veloc_normalized')

flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")


for i in range(num_SNe):
    sampler_8p = sampler.chain[:, :, 8 * i:8 * (i + 1)]
    sampler_1p = sampler.chain[:, :, 8 * num_SNe].reshape(n_walkers,n_steps, 1)
    concat_sampler = np.concatenate((sampler_8p, sampler_1p), axis=2)

    flat_sampler_8p = sampler.get_chain(flat=True)[n_walkers * burn_in:-1, 8 * i:8 * (i + 1)]
    flat_sampler_1p = sampler.get_chain(flat=True)[n_walkers * burn_in:-1, 8 * num_SNe].reshape(n_walkers * (n_steps - burn_in) - 1, 1)
    concat_flat_sampler = np.concatenate((flat_sampler_8p, flat_sampler_1p), axis=1)

    mcmc_snec.plot_lightcurve_with_fit(concat_sampler,
        SN_data_all_w_early[i], parameter_ranges[i],
        'lum_veloc_Tthresh_normalized', res_dir, n_walkers, SN_names[i], n_steps - 1)

    mcmc_snec.plot_lightcurve_with_fit(concat_sampler,
        SN_data_all_w_early[i], parameter_ranges[i],
        'lum_veloc_Tthresh_normalized', res_dir, n_walkers, SN_names[i], 0)
    mcmc_snec.chain_plots(concat_sampler,
                          parameter_ranges[i], res_dir, burn_in, SN_names[i])
    mcmc_snec.corner_plot(concat_flat_sampler,
                          parameter_ranges[i], res_dir, SN_names[i])


