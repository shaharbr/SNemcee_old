import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
# import mcmc_snec
import mcmc_snec_noCSM as mcmc_snec

def loopy_snec_mcmc(SN_name, parameter_ranges):
    extra_tail_days = False  # has to be False or an integer

    Mzams_range = parameter_ranges['Mzams']
    Ni_range = parameter_ranges['Ni']
    E_final_range = parameter_ranges['E']
    # R_range = parameter_ranges['R']
    # K_range = parameter_ranges['K']
    Mix_range = parameter_ranges['Mix']
    T_range = parameter_ranges['T']

    distances = pd.read_csv(os.path.join('results', 'distances.csv'))
    distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
    distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
    sigma_S = 2 * distance_err / distance
    S_range = [1.0 - sigma_S, 1.0 + sigma_S]
    parameter_ranges['S'] = [1.0 - sigma_S, 1.0 + sigma_S]

    # nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

    n_walkers = 200
    n_steps = 1500
    # n_params = 8
    n_params = 6
    burn_in = 500

    time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # res_dir = os.path.join('mcmc_results', str(time_now)+'_lum_'+SN_name)
    res_dir = os.path.join('mcmc_results', str(time_now) + '_lum_noCSM' + SN_name)
    Path(res_dir).mkdir(parents=True, exist_ok=True)

    run_param_df = pd.DataFrame.from_dict({'SN_name': SN_name,
                                 'Mzams_range': str(Mzams_range), 'Ni_range': str(Ni_range),
                                 'E_final_range': str(E_final_range),
                                # 'R_range': str(R_range), 'K_range': str(K_range),
                                 'Mix_range': str(Mix_range), 'Scaling_range': str(S_range), 'T_range': str(T_range),
                                 'n_walkers': n_walkers, 'n_steps': n_steps, 'n_params': n_params,
                                 'burn_in': burn_in,
                                 'time': time_now}, orient='index')
    run_param_df.to_csv(os.path.join(res_dir, 'run_parameters.csv'))

    # import SN bolometric lum data
    data_lum= pd.read_csv(os.path.join('results', SN_name+'_martinez.csv'))

    '''
    # running code #
    '''

    SN_data_all = {'lum': data_lum}

    sampler = mcmc_snec.emcee_fit_params(SN_data_all, n_walkers, n_steps, parameter_ranges,
                                         'lum', extend_tail=extra_tail_days)

    flat_sampler = sampler.get_chain(flat=True)
    np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

    flat_sampler_no_burnin = sampler.get_chain(discard=burn_in, flat=True)
    np.savetxt(os.path.join(res_dir, 'flat_sampler_excluding_burnin.csv'), flat_sampler_no_burnin, delimiter=",")

    mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all, parameter_ranges,
                                            'lum', res_dir, n_walkers, SN_name, n_steps-1, extend_tail=extra_tail_days)
    mcmc_snec.plot_lightcurve_with_fit(sampler.chain, SN_data_all, parameter_ranges,
                                            'lum', res_dir, n_walkers, SN_name, 0, extend_tail=extra_tail_days)

    mcmc_snec.chain_plots(sampler.chain, parameter_ranges, res_dir, burn_in)

    mcmc_snec.corner_plot(sampler.get_chain(flat=True)[n_walkers*burn_in:, :], parameter_ranges, res_dir)

    print(sampler.chain.shape)


Mzams_range = [9.0, 11.0, 13.0, 15.0, 17.0]
Ni_range = [0.02, 0.07, 0.12]
E_final_range = [0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]
# R_range = [0, 500, 1000]
# K_range = [0, 10, 30, 60]
T_range = [-10, 2]
# parameter_ranges_highE = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
#                           'R': R_range, 'K': K_range, 'Mix': Mix_range, 'T': T_range}

parameter_ranges_highE = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                          'Mix': Mix_range, 'T': T_range}

Mzams_range = [9.0, 10.0, 11.0]
Ni_range = [0.001, 0.02, 0.07, 0.12]
E_final_range = [0.1, 0.3, 0.5]
# parameter_ranges_lowE = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
#                           'R': R_range, 'K': K_range, 'Mix': Mix_range, 'T': T_range}

parameter_ranges_lowE = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                          'Mix': Mix_range, 'T': T_range}

SN_list = ['SN2004a', 'SN2004et_d5.9', 'SN2005cs', 'SN2008bk', 'SN2012aw', 'SN2012ec', 'SN2017eaw','SN2018aoq']


for SN_name in SN_list:
    for params in [parameter_ranges_highE]:
        loopy_snec_mcmc(SN_name, params)
