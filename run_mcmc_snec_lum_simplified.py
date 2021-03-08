import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import os
import mcmc_snec_simplified as mcmc_snec

'''
setting the run parameters and loading the luminosity data tables
'''

SN_name = 'SN2012ec'

# MCMC run parameters
n_walkers = 100
n_steps = 2000
burn_in = 500

# user-provided lists of the values for which models exist ("sampled" by computing the model at these values)
Mzams_range = [9.0, 11.0, 13.0, 15.0, 17.0]
Ni_range = [0.02, 0.07, 0.12]
E_final_range = [0.5, 0.9, 1.3, 1.7]
Mix_range = [2.0, 8.0]
R_range = [0, 500, 1000]
K_range = [0, 10, 30, 60]
S_range = [0.5, 2.0]
T_range = [-14, 2]

# organizing the parameter ranges in a dictionary
parameter_ranges = {'Mzams' : Mzams_range, 'Ni' : Ni_range, 'E' : E_final_range,
                    'R' : R_range, 'K' : K_range, 'Mix' : Mix_range, 'S' : S_range, 'T' : T_range}

# defining a non-uniform prior for the S (scaling) parameter, which is set as a gaussian around 1 with uncertainties
# corresponding to the relative uncertainty in distance
distances = pd.read_csv(os.path.join('results', 'distances.csv'))
distance = float(distances.loc[distances['SN_name'] == SN_name]['distance'])
distance_err = float(distances.loc[distances['SN_name'] == SN_name]['distance_err'])
sigma_S = 2 * distance_err / (distance)
nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}

# the name of the folder in which the results of the MCMC will be saved
time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
res_dir = os.path.join('mcmc_results', str(time_now)+'_lum_'+SN_name)
Path(res_dir).mkdir(parents=True, exist_ok=True)

# loading the  SN bolometric lum data
data_lum = pd.read_csv(os.path.join('results', SN_name+'_martinez.csv'))
# take only days between 30 and 196 (like the Martinez paper)
data_lum = data_lum.loc[data_lum['t_from_discovery'] <= 196]


# saving the running parameters and parameter ranges in a csv in the results folder
run_param_df = pd.DataFrame.from_dict({'SN_name': SN_name,
                             'Mzams_range': str(Mzams_range), 'Ni_range': str(Ni_range),
                             'E_final_range': str(E_final_range), 'R_range': str(R_range), 'K_range': str(K_range),
                             'Mix_range': str(Mix_range), 'Scaling_range': str(S_range), 'T_range': str(T_range),
                             'n_walkers': n_walkers, 'n_steps': n_steps,
                             'burn_in': burn_in,
                             'time': time_now}, orient='index')
run_param_df.to_csv(os.path.join(res_dir, 'run_parameters.csv'))


'''
# running the MCMC code #
'''


# run the full MCMC on the lum data and store sampler
sampler = mcmc_snec.emcee_fit_params(data_lum, n_walkers, n_steps, parameter_ranges, nonuniform_priors=nonuniform_priors)

# extract the flat sampler chain (2D matrix, walkersXsteps over parameters) and save as csv
flat_sampler = sampler.get_chain(flat=True)
np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), flat_sampler, delimiter=",")

# plot the final walker's curve fits at the end (last step of the MCMC)
mcmc_snec.plot_lum_with_fit(data_lum, sampler.chain, parameter_ranges, n_steps-1, res_dir, n_walkers, SN_name)
# plot the final walker's curve fits in the beginning (first step of the MCMC)
mcmc_snec.plot_lum_with_fit(data_lum, sampler.chain, parameter_ranges, 0, res_dir, n_walkers, SN_name)
# plot chain plots for all parameters
mcmc_snec.chain_plots(sampler.chain, parameter_ranges, res_dir, burn_in)
# plot corner plot for all parameters, for all walkers over all steps excluding the burn-in
mcmc_snec.corner_plot(sampler.get_chain(flat=True)[n_walkers*burn_in:-1, :], parameter_ranges, res_dir)

print(sampler.chain.shape)
