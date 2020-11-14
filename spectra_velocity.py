import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import random
import data_import
import scipy.integrate as integrate

def approximate_actual_line_range(spectra_df, line_name, lines_dict, absorptionOremission, early=False):
    [range_min, range_max] = lines_dict[line_name][absorptionOremission+'_range']
    if early:
        [range_min, range_max] = [range_min - 100, range_max - 100]
    df_slice = spectra_df.loc[(spectra_df['x'] > range_min) & (spectra_df['x'] < range_max)]
    if absorptionOremission == 'absorption':
        y = np.min(df_slice['y'])
    elif absorptionOremission == 'emission':
        y = np.max(df_slice['y'])
    x = np.min(df_slice.loc[df_slice['y'] == y]['x'])
    width = lines_dict[line_name]['width']
    line_range = [x - width/2, x + width/2]
    return line_range

def find_curve_extreme(polyfit_f, x_slice, param):
    x = list(x_slice)
    y = polyfit_f(x)
    if param == 'absorption':
        extreme_y = np.min(y)
    elif param == 'emission':
        extreme_y = np.max(y)
    extreme_i = np.where(y == extreme_y)[0][0]
    extreme_x = x[extreme_i]
    return {'x': extreme_x, 'y': extreme_y}


def fit_curves(spectra_df, line_name, lines_dict, early, fixed_curve_range=False, number_curves=1):
    # TODO: split this into three seperate functions, for each output
    f, x_slice, curve_extreme = {}, {}, {}
    for param in ['absorption', 'emission']:
        f[param] = list(np.zeros(number_curves))
        x_slice[param] = list(np.zeros(number_curves))
        curve_extreme[param] = list(np.zeros(number_curves))
        if fixed_curve_range:
            line_range = lines_dict[line_name][param+'_range']
        else:
            line_range = approximate_actual_line_range(spectra_df, line_name, lines_dict, param, early=early)
        # calculate multiple polyfits around expected range
        for i in np.arange(0,number_curves):
            range = [random.uniform(line_range[0] - 30, line_range[0] + 30),
                    random.uniform(line_range[1] - 30, line_range[1] + 30)]
            df_slice = spectra_df.loc[(spectra_df['x'] > range[0]) & (spectra_df['x'] < range[1])]
            df_slice.reset_index(inplace=True)
            nans1 = [x for x in df_slice['x'] if any([str(x) == 'nan',  str(x) == 'NaN'])]
            nans2 = [x for x in df_slice['y'] if any([str(x) == 'nan', str(x) == 'NaN'])]
            nans = nans1 + nans2
            if len(df_slice['x']) == 0 or len(nans) > 0:
                f[param][i] = None
                x_slice[param][i] = None
                curve_extreme[param][i] = None
            else:
                slice_poly = np.polyfit(df_slice['x'], df_slice['y'], deg=2)
                f[param][i] = (np.poly1d(slice_poly))
                if (param == 'absorption' and f[param][i][0] > 0) or (param == 'emission' and f[param][i][0] < 0):
                    x_slice[param][i] = (df_slice['x'])
                    curve_extreme[param][i] = find_curve_extreme(f[param][i], x_slice[param][i], param)
                else:
                    f[param][i] = None
                    x_slice[param][i] = None
                    curve_extreme[param][i] = None
    return f, x_slice, curve_extreme



def fit_Pcygni_curves(spectra_dict, lines_dict, fixed_curve_range=False, number_curves=1):
    dates = sorted(list(spectra_dict.keys()))
    lines = lines_dict.keys()
    for date in dates:
        spectra_dict[date]['line'] = {}
        print(spectra_dict[date]['t_from_discovery'])
        if spectra_dict[date]['t_from_discovery'] < 70:
            early = True
        else:
            early = False
        for line_name in lines:
            spectra_dict[date]['line'][line_name] = {}
            spectra_dict[date]['line'][line_name]['polyfits'] , \
            spectra_dict[date]['line'][line_name]['line_xslices'], \
            spectra_dict[date]['line'][line_name]['line_extremes'] \
                = fit_curves(spectra_dict[date]['df'], line_name, lines_dict, early, fixed_curve_range, number_curves)
    return spectra_dict



def add_fitted_curves_to_plot(single_spectrum_dict, lines_dict, ax, y_shift, line_velocity, number_curves_to_draw=1):
    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#d95f02', 'FeII 5018': '#1D78EF'}
    if line_velocity:
        lines = [line_velocity]
    else:
        lines = lines_dict.keys()
    for line_name in lines:
        line_dict = single_spectrum_dict['line'][line_name]
        for i in range(number_curves_to_draw):
            f = line_dict['polyfits']['absorption'][i]
            if f != None:
                x_slice = line_dict['line_xslices']['absorption'][i]
                curve_y = f(x_slice) + y_shift
                extreme_point = line_dict['line_extremes']['absorption'][i]
                x_point = extreme_point['x']
                y_point = extreme_point['y'] + y_shift
                if line_velocity:
                    wavelength_expected = lines_dict[line_name]['peak']
                    x_slice = calculate_expansion_velocity(wavelength_expected, x_slice)
                    x_point = calculate_expansion_velocity(wavelength_expected, [x_point])
                ax.plot(x_slice, curve_y, color=colors[line_name], alpha=0.1, linewidth=1.5)
                ax.plot(x_point, y_point, '.', color=colors[line_name])


def calculate_expansion_velocity(wavelength_expected, wavelength_observed):
    num_curves = len(wavelength_observed)
    wavelength_expected = [wavelength_expected] * num_curves
    c = 299792.458  # km/s
    v = c * (np.divide(wavelength_expected, wavelength_observed) - 1)
    return v

def add_expansion_velocity(spectra_dict, lines_dict):
    dates = sorted(list(spectra_dict.keys()))
    lines = lines_dict.keys()
    wavelength_observed = {}
    for date in dates:
        for line_name in lines:
            # because of the blue-ness of the early spectra, any Pcygni other than Halpha
            # are very hard to identify (and therefore innacurate)
            if line_name == 'Halpha' or spectra_dict[date]['t_from_discovery'] > 20:
                spectra_dict[date]['line'][line_name]['velocity'] = {}
                wavelength_expected = lines_dict[line_name]['peak']
                for param in ['absorption', 'emission']:
                    spectra_dict[date]['line'][line_name]['velocity'][param] = {}
                    number_curves = len(spectra_dict[date]['line'][line_name]['line_extremes'][param])
                    print(spectra_dict[date]['line'][line_name]['line_extremes'][param])
                    # if None not in spectra_dict[date]['line'][line_name]['line_extremes'][param]:
                    wavelength_observed[param] = []
                    for i in range(number_curves):
                        if spectra_dict[date]['line'][line_name]['line_extremes'][param][i] != None:
                            wavelength_observed[param].append(
                                spectra_dict[date]['line'][line_name]['line_extremes'][param][i]['x'])
                    velocity_list = calculate_expansion_velocity(wavelength_expected, wavelength_observed[param])
                    spectra_dict[date]['line'][line_name]['velocity'][param]['mean'] = np.mean(velocity_list)
                    print(spectra_dict[date]['t_from_discovery'])
                    print(line_name)
                    print(param)
                    print(velocity_list)
                    print(np.std(velocity_list))
                    spectra_dict[date]['line'][line_name]['velocity'][param]['std'] = np.std(velocity_list)
                    if len(wavelength_observed[param]) < 1:
                        print('none active')
                        spectra_dict[date]['line'][line_name]['velocity'][param]['mean'] = None
                        spectra_dict[date]['line'][line_name]['velocity'][param]['std'] = None
    return spectra_dict

def add_name_t_from_discovery_to_df(df, name, discovery_date):
    df['datetime'] = data_import.convert_to_mjd(df['datetime'], from_datetime=True)
    df['t_from_discovery'] = df['datetime'] - discovery_date
    df['Name'] = name
    return df


def make_velocity_df(SN_dict, lines_dict):
    dates = SN_dict['spectra'].keys()
    lines = lines_dict.keys()
    velocity_df = []
    for date in dates:
        for line in lines:
            # because of the blue-ness of the early spectra, any Pcygni other than Halpha
            # are very hard to identify (and therefore innacurate)
            if line == 'Halpha' or SN_dict['spectra'][date]['t_from_discovery'] > 20:
                velocity_df.append(pd.DataFrame({'line': line,
                                                 'absorption_mean_velocity': SN_dict['spectra'][date]['line'][line]['velocity']['absorption']['mean'],
                                                 'absorption_std_velocity': SN_dict['spectra'][date]['line'][line]['velocity']['absorption']['std'],
                                                 'emission_mean_velocity': SN_dict['spectra'][date]['line'][line]['velocity']['emission']['mean'],
                                                 'emission_std_velocity':SN_dict['spectra'][date]['line'][line]['velocity']['emission']['std'],
                                                 't_from_discovery': SN_dict['spectra'][date]['t_from_discovery']},
                                                  index=[0]))
    velocity_df = pd.concat(velocity_df)
    velocity_df['Name'] = SN_dict['Name']
    return velocity_df


# plot expansion velocities over datetime
def plot_expansion_velocities(df_list, absorptionORemission):
    fig, ax = plt.subplots(1, figsize=(9, 6))

    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#d95f02', 'FeII 5018': '#1D78EF'}
    markers_fill = {'SN2018hmx': 'full', 'SNiPTF14hls': 'none', 'SN2018aad': 'none'}
    names = []
    for df in df_list:
        SN_name = pd.unique(df['Name'])[0]
        names.append(SN_name)
        lines = colors.keys()
        for line_name in lines:
            linedf = df.loc[df['line'] == line_name]
            ax.errorbar(x=linedf['t_from_discovery'], y=linedf[absorptionORemission+'_mean_velocity'],
                          yerr=linedf[absorptionORemission+'_std_velocity'],
                                                 label=SN_name+' '+line_name,
                                                 marker='s',
                                                 fillstyle=markers_fill[SN_name],
                                                 linestyle='None',
                                                 color=colors[line_name])




    if absorptionORemission == 'absorption':
        ax.set_ylim(bottom=0)
    ax.set_title('Expansion velocity over time - '+absorptionORemission + ' - '+ str(names))
    ax.set_ylabel('Expansion velocity (km/s)')
    ax.set_xlabel('Rest-frame days from discovery')
    ax.legend()
    ax.tick_params(axis='both', which='major')
    fig.savefig(r'figures/' + ''.join(names) + absorptionORemission + '_expansion_velocity_over_time_sn1999em' + '.png')
    fig.savefig(r'figures/' + ''.join(names) + absorptionORemission + '_expansion_velocity_over_time_sn1999em' + '.svg')



def pEW(spectra_dict, line_name, date):
    result = []
    # plt.figure()
    t_from_discovery = spectra_dict[date]['t_from_discovery']
    # plt.title('day '+ str(int(t_from_discovery)))
    for i in np.arange(20):
        # print(spectra_dict[date])
        # print(spectra_dict.keys())
        # print(spectra_dict[date].keys())
        # print(spectra_dict[date].keys())
        line_xslice = spectra_dict[date]['line'][line_name]['line_xslices']['absorption'][i]
        if line_xslice is not None:
            f_polyfit = spectra_dict[date]['line'][line_name]['polyfits']['absorption'][i]
            line_xslice = list(line_xslice)
            spectra_df = spectra_dict[date]['df']
            fit_slice = []
            x_left = line_xslice[0]
            x_right = line_xslice[-1]
            for x in [x_left, x_right]:
                x_min, x_max = x - 20, x + 20
                df_slice = spectra_df.loc[(spectra_df['x'] > x_min) & (spectra_df['x'] < x_max)]
                y_peak = np.max(df_slice['y'])
                x = np.min(df_slice.loc[df_slice['y'] == y_peak]['x'])
                x_min, x_max = x - 5, x + 5
                df_slice = spectra_df.loc[(spectra_df['x'] > x_min) & (spectra_df['x'] < x_max)]
                fit_slice.append(df_slice)
                range_scatter = spectra_df.loc[(spectra_df['x'] > x_min) & (spectra_df['x'] < x_max)]
                # plt.scatter(range_scatter['x'], range_scatter['y'], color='blue', alpha=0.1)
            fit_slice = pd.concat(fit_slice, ignore_index=True)
            line_fit = np.poly1d(np.polyfit(fit_slice['x'], fit_slice['y'], deg=1))

            range_scatter = spectra_df.loc[(spectra_df['x'] > x_left - 50) & (spectra_df['x'] < x_right + 50)]
            # plt.scatter(range_scatter['x'], range_scatter['y'], color='orange', alpha=0.1)

            common_x = list(range_scatter['x'])
            line_y = line_fit(common_x)
            curve_y = f_polyfit(common_x)
            # plt.plot(common_x, line_y, 'k', alpha=0.1)
            # plt.plot(common_x, curve_y, '#1D78EF', alpha=0.1)


            intersect_i = np.argwhere(np.diff(np.sign(line_y - curve_y))).flatten()
            x_intersect = [common_x[i] for i in intersect_i]
            y_intersect = [line_fit(x) for x in x_intersect]
            # plt.plot(x_intersect, y_intersect, 'ro', alpha=0.1)
            # plt.xlabel('Rest (z = 0.038) Wavelength (Å)')
            # plt.ylabel('Normalized fλ')
            if len(x_intersect) > 1:
                result.append(integrate.quad(lambda x: 1 - (f_polyfit(x) / (line_fit(x))), x_intersect[0], x_intersect[1])[0])
    n = len(result)
    result_mean = np.mean(result)
    result_std = np.std(result)
    # print('final', str(t_from_discovery), 'days',  result_mean, result_std, 'n=', n)
    return result_mean, result_std
        # plt.plot(common_x, curve_y /line_y)
        # plt.plot(common_x, 1 - (curve_y / (line_y + 0.1)))
        # div = curve_y /line_y
        # pEW_func = f_polyfit / line_fit
        # print('f_polyfit', f_polyfit)
        # print('line_fit', line_fit)
        # print('pEW_func', pEW_func)
        # plt.plot(common_x, pEW_func(common_x))



        # number = 0
        # return number
