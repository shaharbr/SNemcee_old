import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import mcmc_snec
import re
import snec_model_interpolator as interp
import corner

T_thresh = 10 ** 3.75
result_name = '2021-04-16_16-37-35_lum_veloc_SN2004a'
result_dir = os.path.join('mcmc_results', result_name)
extend_tail = False
models = {}


def initialize_empty_models_dict(ranges_dict):
    for data_type in ['lum', 'veloc', 'mag', 'temp']:
        models[data_type] = {}
        for Mzams in ranges_dict['Mzams']:
            models[data_type][Mzams] = {}
            for Ni in ranges_dict['Ni']:
                models[data_type][Mzams][Ni] = {}
                for E in ranges_dict['E']:
                    models[data_type][Mzams][Ni][E] = {}
                    for R in ranges_dict['R']:
                        models[data_type][Mzams][Ni][E][R] = {}
                        for K in ranges_dict['K']:
                            models[data_type][Mzams][Ni][E][R][K] = {}
                            for Mix in ranges_dict['Mix']:
                                models[data_type][Mzams][Ni][E][R][K][Mix] = None


def get_surrouding_values(requested, ranges_dict):
    params = list(ranges_dict.keys())
    surrouding_values = {param: [] for param in params}
    for i in range(len(requested)):
        param_range = np.array(ranges_dict[params[i]])
        below = np.max(param_range[param_range <= requested[i]])
        above = np.min(param_range[param_range >= requested[i]])
        surrouding_values[params[i]] = [below, above]
    return surrouding_values


def load_model(Mzams, Ni, E, R, K, Mix, data_type, extend_tail=False):
    if R == 0 or K == 0:
        R = 0
        K = 0
    name = 'M' + str(Mzams) + \
           '_Ni' + str(Ni) + \
           '_E' + str(E) + \
           '_Mix' + str(Mix) + \
           '_R' + str(R) + \
           '_K' + str(K)
    if data_type == 'lum':
        modelpath = os.path.join('..', 'all_lum_data', name, 'lum_observed.dat')
        if os.stat(modelpath).st_size < 10 ** 5:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath,
                                     names=['t_from_discovery', 'Lum'], sep=r'\s+')
            time_col = snec_model['t_from_discovery'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, time_col, snec_model['Lum'])
            if extend_tail is not False:
                last_30d_x = interp_days[-100:-1]
                last_30d_y = snec_model[-100:-1]
                last_30d_ylog = np.log(last_30d_y)
                tail_poly1d = np.poly1d(np.polyfit(last_30d_x, last_30d_ylog, deg=1))
                extension_days = np.linspace(200.1, 200+extend_tail, int(10*extend_tail))
                extension_lumlog = np.array([tail_poly1d(extension_days[i]) for i in range(len(extension_days))])
                extension_lum = np.exp(extension_lumlog)
                snec_model = np.concatenate((snec_model, extension_lum))
            return snec_model
    elif data_type == 'veloc':
        modelpath = os.path.join('..', 'all_veloc_data', name, 'vel_Fe.dat')
        if os.stat(modelpath).st_size < 10 ** 4:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath)
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, snec_model['t_from_discovery'],
                                   snec_model['veloc'])
            return snec_model
    elif data_type == 'temp':
        modelpath = os.path.join('..', 'all_temp_rad_data', name, 'T_eff.dat')
        if os.stat(modelpath).st_size < 10 ** 5:
            return 'failed SN'
        snec_model = pd.read_csv(modelpath,
                                 names=['t_from_discovery', 'temp'], sep=r'\s+')
        time_col = snec_model['t_from_discovery'] / 86400  # sec to days
        interp_days = np.linspace(0, 200, 2001)
        snec_model = np.interp(interp_days, time_col, snec_model['temp'])
        return snec_model
    elif data_type == 'mag':
        modelpath = os.path.join('..', 'all_pys_mag_data', name, 'magnitudes.dat')
        lumpath = os.path.join('..', 'all_lum_data', name, 'lum_observed.dat')
        if os.stat(lumpath).st_size < 10 ** 5:
            return 'failed SN'
        else:
            mag_file = pd.read_csv(modelpath,
                                   names=['time', 'Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'],
                                   sep=r'\s+')
            mag_file = mag_file.abs()
            time_col = mag_file['time'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model_dict = {}
            for filter in ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
                snec_model_dict[filter] = np.interp(interp_days, time_col, mag_file[filter])
            snec_model_dict['time'] = interp_days
            snec_model = pd.DataFrame(snec_model_dict)
            snec_model = snec_model.sort_values('time')
            return snec_model


def load_surrounding_models(requested, ranges_dict, fitting_type, extend_tail=False):
    surrouding_values = get_surrouding_values(requested, ranges_dict)
    for Mzams in surrouding_values['Mzams']:
        for Ni in surrouding_values['Ni']:
            for E in surrouding_values['E']:
                for R in surrouding_values['R']:
                    for K in surrouding_values['K']:
                        for Mix in surrouding_values['Mix']:
                            if 'lum' in fitting_type:
                                if models['lum'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['lum'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'lum', extend_tail)
                            if 'veloc' in fitting_type:
                                if models['veloc'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['veloc'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'veloc')
                                if models['temp'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['temp'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'temp')
                            if 'mag' in fitting_type:
                                if models['mag'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['mag'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'mag')



def import_ranges(list):
    ranges = []
    for name in list:
        if 'range' in name:
            content = re.sub(r'[\[\]\s]', '', run_params.iloc[0][name]).split(',')
            if name == 'R_range' or name == 'K_range':
                content = [int(content[i]) for i in range(len(content))]
            else:
                content = [float(content[i]) for i in range(len(content))]
            ranges.append(content)
    return ranges


def get_param_results_dict(sampler_df):
    dict = {}
    for param in ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']:
        avg = np.average(sampler_df[param])
        sigma_lower, sigma_upper = np.percentile(sampler_df[param], [16, 84])
        dict[param] = avg
        dict[param + '_lower'] = avg - sigma_lower
        dict[param + '_upper'] = sigma_upper - avg
    return dict


def result_text_from_dict(sampler_df):
    param_dict = get_param_results_dict(sampler_df)
    res_text = 'MCMC-SNEC fit\n\n'
    for param in ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']:
        if (param != 'K') & (param != 'R'):
            res_text += param + ': ' + mcmc_snec.rounded_str(param_dict[param]) + r'$\pm$ [' +\
                            mcmc_snec.rounded_str(param_dict[param+'_lower']) + ',' +\
                            mcmc_snec.rounded_str(param_dict[param+'_upper']) + ']\n'
    if (param_dict['K'] == 0) & (param_dict['R'] == 0):
        res_text += 'no CSM'
    else:
        res_text += 'K' + ': ' + mcmc_snec.rounded_str(param_dict['K']) + r'$\pm$ [' + \
                    mcmc_snec.rounded_str(param_dict['K_lower']) + ',' + \
                    mcmc_snec.rounded_str(param_dict['K_upper']) + ']\n'
        res_text += 'R' + ': ' + mcmc_snec.rounded_str(param_dict['R']) + r'$\pm$ [' + \
                    mcmc_snec.rounded_str(param_dict['R_lower']) + ',' + \
                    mcmc_snec.rounded_str(param_dict['R_upper']) + ']\n'
    return res_text


def plot_lum_with_fit(data_lum, sampler_df, ranges_dict, output_dir, n_walkers, SN_name, martinez_values, extend_tail=False):
    data_x = data_lum['t_from_discovery']
    data_y = data_lum['Lum']
    dy0 = data_lum['dLum0']
    dy1 = data_lum['dLum1']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    add_label = True
    if extend_tail is not False:
        x_plotting = np.linspace(0, 200 + extend_tail, int(1 + 10 * 200 + extend_tail))
    else:
        x_plotting = np.linspace(0, 200, 2001)
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict) or R == 0:
            data_x_moved = x_plotting - T
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            load_surrounding_models(requested[0:6], ranges_dict, 'lum', extend_tail)
            y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['lum'], data_x_moved, extend_tail)
            if not isinstance(y_fit, str):
                # multiply whole graph by scaling factor
                y_fit = y_fit * S
                if add_label:
                    ax.plot(x_plotting, np.log10(y_fit), alpha=0.05, color='purple', label='SNEC-MCMC fits')
                    add_label = False
                else:
                    ax.plot(x_plotting, np.log10(y_fit), alpha=0.05, color='purple')
    data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')

    results_text = result_text_from_dict(sampler_df)
    ax.text(1.02, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='purple')
    # martinez fits (according to SNEC)
    data_x_moved = x_plotting - martinez_values[7]
    surrounding_values = get_surrouding_values(martinez_values[0:6], ranges_dict)
    load_surrounding_models(martinez_values[0:6], ranges_dict, 'lum', extend_tail)
    y_fit = interp.snec_interpolator(martinez_values[0:6], surrounding_values, models['lum'], data_x_moved, extend_tail)
    y_fit = y_fit * martinez_values[6]
    ax.plot(x_plotting, np.log10(y_fit), color='green', label='Martinez values on SNEC', linewidth=2.0)

    # martinez fits (according to their model)
    data_x_moved = martinez_fits['TIME'] + martinez_values[7]
    y_fit = np.power(10, martinez_fits['LOGL'])
    y_fit = y_fit * martinez_values[6]
    ax.plot(data_x_moved, np.log10(y_fit), color='orange', label='Martinez fits', linewidth=2.0)

    results_text = 'Martinez fits' + '\n\n' + \
                   'Mzams: ' + str(martinez_values[0]) + '\n' + \
                   'Ni: ' + str(martinez_values[1]) + '\n' + \
                   'E: ' + str(martinez_values[2]) + '\n' + \
                   'Mix: ' + str(martinez_values[5]) + '\n' + \
                   'S: ' + str(martinez_values[6]) + '\n' + \
                   'T: ' + str(martinez_values[7]) + '\n' + \
                   'noCSM'
    ax.text(1.02, 0.4, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='green')
    # real observations for the SN
    ax.errorbar(data_x, data_y, yerr=[dy0, dy1], marker='o', linestyle='None', color='k',
                label=SN_name + ' observations')
    ax.set_xlim(-2, 200)
    ax.set_title(SN_name, fontsize=18)
    ax.legend()
    plt.tight_layout()
    if extend_tail is not False:
        ax.set_xlim(-2, 200 + extend_tail)
    else:
        ax.set_xlim(-2, 200)
    ax.set_ylim(40.8, 43)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, str(result_name) + '_lum.png'))


def plot_veloc_with_fit(data_veloc, sampler_df, ranges_dict, output_dir, n_walkers, SN_name, martinez_values, Tthreshold=False):
    data_x = data_veloc['t_from_discovery']
    data_y = data_veloc['veloc']
    data_dy = data_veloc['dveloc']
    f_fit, ax = plt.subplots(figsize=(10, 8))
    x_plotting = np.linspace(0, np.max(data_veloc['t_from_discovery']), int(1+10*np.max(data_veloc['t_from_discovery'])))
    add_label = True
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict) or R == 0:
            surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
            load_surrounding_models(requested[0:6], ranges_dict, 'veloc')
            temp_x_moved = data_x - T
            if Tthreshold:
                temp_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['temp'], temp_x_moved)
                if not isinstance(temp_fit, str):
                    max_temp_below_Tthresh = np.max(temp_fit[temp_fit <= T_thresh])
                    temp_x_moved = temp_x_moved[temp_fit > max_temp_below_Tthresh]
                    if len(temp_x_moved) <= 1:
                        print('cooled too fast, no early velocity data')
                        x_plotting = []
                    else:
                        max_x_moved = np.max(temp_x_moved)
                        x_plotting = np.linspace(0, max_x_moved, int(1+max_x_moved*10))
            else:
                x_plotting = np.linspace(0, 200, 2001)
            if len(x_plotting) > 0:
                data_x_moved = x_plotting - T
                y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['veloc'], data_x_moved)
                if not isinstance(y_fit, str):
                    if add_label:
                        ax.plot(x_plotting, y_fit, alpha=0.05, color='purple', label='SNEC-MCMC fits')
                        add_label = False
                    else:
                        ax.plot(x_plotting, y_fit, alpha=0.05, color='purple')

    results_text = result_text_from_dict(sampler_df)
    ax.text(1.02, 0.8, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='purple')
    # martinez fits (according to SNEC)
    data_x_moved = x_plotting - martinez_values[7]
    surrounding_values = get_surrouding_values(martinez_values[0:6], ranges_dict)
    load_surrounding_models(martinez_values[0:6], ranges_dict, 'veloc')
    y_fit = interp.snec_interpolator(martinez_values[0:6], surrounding_values, models['veloc'], data_x_moved)
    ax.plot(x_plotting, y_fit, color='green', label='Martinez values on SNEC', linewidth=2.0)

    # martinez fits (according to their model)
    data_x_moved = martinez_fits['TIME'] + martinez_values[7]
    ax.plot(data_x_moved, martinez_fits['VPH'], color='orange', label='Martinez fits', linewidth=2.0)

    results_text = 'Martinez fits' + '\n\n' + \
                   'Mzams: ' + str(martinez_values[0]) + '\n' + \
                   'Ni: ' + str(martinez_values[1]) + '\n' + \
                   'E: ' + str(martinez_values[2]) + '\n' + \
                   'Mix: ' + str(martinez_values[5]) + '\n' + \
                   'S: ' + str(martinez_values[6]) + '\n' + \
                   'T: ' + str(martinez_values[7]) + '\n' + \
                   'noCSM'
    ax.text(1.02, 0.4, results_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='center', bbox=dict(facecolor='white', alpha=0.5), color='green')
    # real observations for the SN
    ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', color='k',
                label=SN_name + ' observations')
    ax.legend()
    ax.set_xlim(-2, 140)
    ax.set_title(SN_name, fontsize=18)
    plt.tight_layout()
    f_fit.savefig(os.path.join(output_dir, str(result_name) + '_veloc.png'))


def range_bounds(ranges_list):
    tup_list = []
    for i in range(len(ranges_list)):
        tup_list.append((np.min(ranges_list[i]), np.max(ranges_list[i])))
    return tup_list


def corner_plot(sampler_chain_flat, ranges_list, output_dir, martinez_values):
    labels = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']
    num_param = len(labels)
    corner_range = range_bounds(ranges_list)
    f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
    # Extract the axes
    axes = np.array(f_corner.axes).reshape((num_param, num_param))
    # Loop over the diagonal
    for i in range(num_param):
        ax = axes[i, i]
        ax.axvline(martinez_values[i], color="orange")
    # Loop over the histograms
    for yi in range(num_param):
        for xi in range(yi):
            ax = axes[yi, xi]
            ax.axvline(martinez_values[xi], color="orange")
            ax.axhline(martinez_values[yi], color="orange")
            ax.plot(martinez_values[xi], martinez_values[yi], marker='s', markersize=10, color='darkorange')
    f_corner.savefig(os.path.join(output_dir, str(result_name) + '_corner_plot.png'))

flat_sampler_path = os.path.join(result_dir, 'flat_sampler.csv')
output_path = os.path.join(result_dir, 'run_parameters.csv')

run_params_path = os.path.join(result_dir, 'run_parameters.csv')
run_params = pd.read_csv(run_params_path, index_col=0).T

ranges = import_ranges(run_params.columns.values)
ranges_dict = {'Mzams' : ranges[0], 'Ni' : ranges[1], 'E' : ranges[2],
               'R' : ranges[3], 'K' : ranges[4], 'Mix' : ranges[5],
               'S' : ranges[6], 'T' : ranges[7]}

SN_name = run_params.iloc[0]['SN_name']
n_walkers = int(run_params.iloc[0]['n_walkers'])
n_steps = int(run_params.iloc[0]['n_steps'])
burn_in = int(run_params.iloc[0]['burn_in'])



flat_sampler_all = pd.read_csv(flat_sampler_path,
                               names=['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T'],
                               skiprows=(burn_in-1) * (n_walkers))

flat_sampler = pd.read_csv(flat_sampler_path,
                               names=['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T'],
                               skiprows=(n_steps-1) * (n_walkers))


martinez_values_path = os.path.join('results', SN_name+'_martinez_values.csv')
martinez_values = list(pd.read_csv(martinez_values_path, names=['val']).T.iloc[0].astype(float))

martinez_fits_path = os.path.join('results', 'martinez_bestfit_models', SN_name+'.dat')
martinez_fits = pd.read_csv(martinez_fits_path, sep=r'\s+')
martinez_fits['VPH'] = martinez_fits['VPH'] / 10**5   # cm/s to km/s

data_lum_path = os.path.join('results', SN_name+'_martinez.csv')
data_veloc_path = os.path.join('results', SN_name+'_expansion_velocities.csv')
data_lum = pd.read_csv(data_lum_path)
data_veloc = pd.read_csv(data_veloc_path,
                      usecols=['t_from_discovery', 'line', 'absorption_mean_velocity','absorption_std_velocity'])
data_veloc.rename({'absorption_mean_velocity':'veloc', 'absorption_std_velocity':'dveloc'}, axis='columns', inplace=True)
data_veloc = data_veloc.loc[data_veloc['line'] == 'FeII 5169']
data_veloc.sort_values(by=['t_from_discovery'], inplace=True)
# remove first point which seems like an artifact
data_veloc = data_veloc.loc[data_veloc['t_from_discovery'] > 20]

SN_data = {'lum': data_lum, 'veloc': data_veloc}

initialize_empty_models_dict(ranges_dict)

# #### running functions ################################

corner_plot(flat_sampler_all, ranges, 'figures', martinez_values)

plot_lum_with_fit(SN_data['lum'], flat_sampler, ranges_dict, 'figures', n_walkers, SN_name, martinez_values, extend_tail)

plot_veloc_with_fit(SN_data['veloc'], flat_sampler, ranges_dict, 'figures', n_walkers, SN_name, martinez_values, T_thresh)


