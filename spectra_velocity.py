import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import re
import random
import data_import

plt.rcParams['font.sans-serif'] = 'Arial'





def add_rest_frame_days_from_discovery(SN_dict):
    dates = SN_dict['spectra'].keys()
    z = SN_dict['z']
    for date in dates:
        SN_dict['spectra'][date]['t_from_discovery'] = (date - SN_dict['discovery_date']) / (1 + z)
    return SN_dict


def correct_specta_redshift(SN_dict):
    dates = SN_dict['spectra'].keys()
    z = SN_dict['z']
    for date in dates:
        # correct wavelengths for redshift
        SN_dict['spectra'][date]['df']['x'] = SN_dict['spectra'][date]['df']['x'] / (1 + z)
    return SN_dict

def normalize_spectra(SN_dict):
    dates = SN_dict['spectra'].keys()
    for date in dates:
        avg = np.average(SN_dict['spectra'][date]['df']['y'])
        SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / avg
        if SN_dict['spectra'][date]['telescope'] == 'ZTF':
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] * 2
        if SN_dict['spectra'][date]['telescope'] == 'HET':
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / 2
    return SN_dict



def smooth_spectrum(spectrum_y):
    var_list = list(spectrum_y)
    var_list = [x for x in var_list if str(x) != 'nan']
    var_list = var_list[0:200]
    start_variability = np.std(var_list)
    smoothing_window = int(len(spectrum_y)/50 * start_variability)+1
    smooth_y = spectrum_y.rolling(smoothing_window, center=True).mean()
    return smooth_y


# plot all spectra in an overlay figure
def plot_overlay_spectra(SN_dict, lines_dict=False):
    dates = SN_dict['spectra'].keys()
    fig, ax = plt.subplots(1, figsize=(10, 5))
    name = SN_dict['name']
    for date in dates:
        SN_dict['spectra'][date]['df'].plot(x='x', y='y', ax=ax, label=date)
    if lines_dict:
        for line_name in lines_dict.keys():
            ax.axvline(x=lines_dict[line_name]['peak'], color='k')
    ax.set_title(name+' spectra over time', fontsize=16)
    ax.set_xlabel('Rest (z = 0.037) Wavelength (Å)', size=14)
    ax.set_ylabel('Flux (10-15 erg s-1 cm-2 Å-1)', size=14)
    ax.tick_params(axis='both', which='major', labelsize=14)
    # fig.savefig(name + '_spectra_over_time_overlay'+'.png')


def approximate_actual_line_range(spectra_df, line_name, lines_dict, absorptionOremission):
    [range_min, range_max] = lines_dict[line_name][absorptionOremission+'_range']
    df_slice = spectra_df.loc[(spectra_df['x'] > range_min) & (spectra_df['x'] < range_max)]
    if absorptionOremission == 'absorption':
        y = np.min(df_slice['y'])
    elif absorptionOremission == 'emission':
        y = np.max(df_slice['y'])
    x = np.min(df_slice.loc[df_slice['y'] == y]['x'])
    line_range = [x - 80, x + 80]
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


def fit_curves(spectra_df, line_name, lines_dict, fixed_curve_range=False, number_curves=1):
    # TODO: split this into three seperate functions, for each output
    f, x_slice, curve_extreme = {}, {}, {}
    for param in ['absorption', 'emission']:
        f[param] = list(np.zeros(number_curves))
        x_slice[param] = list(np.zeros(number_curves))
        curve_extreme[param] = list(np.zeros(number_curves))
        if fixed_curve_range:
            line_range = lines_dict[line_name][param+'_range']
        else:
            line_range = approximate_actual_line_range(spectra_df, line_name, lines_dict, param)
        # calculate multiple polyfits around expected range
        for i in np.arange(0,number_curves):
            range = [random.uniform(line_range[0] - 5, line_range[0] + 5),
                    random.uniform(line_range[1] - 5, line_range[1] + 5)]
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
    dates = spectra_dict.keys()
    lines = lines_dict.keys()
    for date in dates:
        spectra_dict[date]['line'] = {}
        for line_name in lines:
            spectra_dict[date]['line'][line_name] = {}
            spectra_dict[date]['line'][line_name]['polyfits'] , \
            spectra_dict[date]['line'][line_name]['line_xslices'], \
            spectra_dict[date]['line'][line_name]['line_extremes']\
                = fit_curves(spectra_dict[date]['df'], line_name, lines_dict, fixed_curve_range, number_curves)
    return spectra_dict


def plot_stacked_spectra(SN_dict, lines_dict, plot_curve_fits=False, line_velocity=False):
    dates = sorted(SN_dict['spectra'].keys(), reverse=True)
    if line_velocity:
        fig, ax = plt.subplots(1, figsize=(5, 9))
    else:
        fig, ax = plt.subplots(1, figsize=(12, 9))
    name = SN_dict['name']
    y_shift = 0
    lines_list = lines_dict.keys()
    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#''d95f02'}
    # draw expected lines locations in background
    if line_velocity:
        ax.axvline(x=0, color=colors[line_velocity], alpha=0.5)
    else:
        for line_name in lines_list:
            wavelength_expected = lines_dict[line_name]['peak']
            ax.axvline(x=wavelength_expected, color=colors[line_name], alpha=0.5)
    # draw spectra
    for date in dates:
        date_dict = SN_dict['spectra'][date]
        df = date_dict['df']
        if date_dict['t_from_discovery'] < 20:
            y_shift = y_shift - 0.6
            smooth_y = df['y']
        else:
            smooth_y = smooth_spectrum(df['y'])
        if line_velocity:
            if line_velocity != 'Halpha':
                y_shift_delta = 0.3
                if date_dict['t_from_discovery'] < 20:
                    y_shift = y_shift - 0.6
            else:
                y_shift_delta = 0.8
            min_x = lines_dict[line_velocity]['absorption_range'][0] - 150
            max_x = lines_dict[line_velocity]['peak'] + 200
            bool_x = (df['x'] > min_x) & (df['x'] < max_x)
            x = df['x'].loc[bool_x]
            y = smooth_y[bool_x]
            # y = df['y'].loc[bool_x]
            wavelength_expected = lines_dict[line_velocity]['peak']
            x = calculate_expansion_velocity(wavelength_expected, x)
            max_x = np.nanmax(x)
            min_x = np.nanmin(x)
            y = y + y_shift
            ax.set_xlim(min_x, max_x)
            label_x = min_x - 150
            ax.invert_xaxis()
            ax.set_title(line_velocity, fontsize=16)
            ax.set_xlabel('Velocity (km/s)', size=14)
        else:
            x = df['x']
            y_shift_delta = (np.max(smooth_y) - np.min(smooth_y)) / 3
            y = smooth_y + y_shift
            ax.set_xlim(3300, 9700)
            label_x = 9780
            ax.set_title(name + ' spectra over time', fontsize=16)
            ax.set_xlabel('Rest (z = 0.037) Wavelength (Å)', size=14)
        # plot
        ax.plot(x, y, color='navy', linewidth=1)
        # spectrum labeled with day from explosion
        spectrum_date_label = str(int(date_dict['t_from_discovery'])) + 'd (' + date_dict['telescope'] + ')'
        label_ypos = np.mean(y) - 0.4
        ax.text(label_x, label_ypos, spectrum_date_label , fontsize=14, color='navy', fontname="Arial")
        # adding the polyfit curves if asked
        if plot_curve_fits:
            add_fitted_curves_to_plot(date_dict, lines_dict, ax, y_shift, line_velocity, number_curves_to_draw=3)
        y_shift += y_shift_delta
    ax.set_ylabel('Normalized fλ + shift', size=14)
    ax.set_yticks([])
    # ax.get_legend().remove()
    ax.tick_params(axis='both', which='major', labelsize=14)
    fig.subplots_adjust(right=0.7)
    if plot_curve_fits:
        fig.savefig(name + '_spectra_over_time_stacked_with_polyfit' + str(line_velocity) + '.png')
    else:
        fig.savefig(name + '_spectra_over_time_stacked' + str(line_velocity) + '.png')



def add_fitted_curves_to_plot(single_spectrum_dict, lines_dict, ax, y_shift, line_velocity, number_curves_to_draw=1):
    colors = {'Halpha':'r', 'Hbeta':'g', 'FeII 5169':'b'}
    if line_velocity:
        lines = [line_velocity]
    else:
        lines = lines_dict.keys()
    for line_name in lines:
        line_dict = single_spectrum_dict['line'][line_name]
        for param in ['absorption']:
            for i in range(number_curves_to_draw):
                f = line_dict['polyfits'][param][i]
                if f != None:
                    x_slice = line_dict['line_xslices'][param][i]
                    curve_y = f(x_slice) + y_shift
                    extreme_point = line_dict['line_extremes'][param][i]
                    x_point = extreme_point['x']
                    y_point = extreme_point['y'] + y_shift
                    if line_velocity:
                        wavelength_expected = lines_dict[line_name]['peak']
                        x_slice = calculate_expansion_velocity(wavelength_expected, x_slice)
                        x_point = calculate_expansion_velocity(wavelength_expected, [x_point])
                    ax.plot(x_slice, curve_y, color=colors[line_name], alpha=0.5, linewidth=1.5)
                    ax.plot(x_point, y_point, '.', color=colors[line_name])


def calculate_expansion_velocity(wavelength_expected, wavelength_observed):
    num_curves = len(wavelength_observed)
    wavelength_expected = [wavelength_expected] * num_curves
    c = 299792.458  # km/s
    v = c * (np.divide(wavelength_expected, wavelength_observed) - 1)
    return v

def add_expansion_velocity(spectra_dict, lines_dict):
    dates = spectra_dict.keys()
    lines = lines_dict.keys()
    wavelength_observed = {}
    for date in dates:
        for line_name in lines:
            spectra_dict[date]['line'][line_name]['velocity'] = {}
            wavelength_expected = lines_dict[line_name]['peak']
            for param in ['absorption', 'emission']:
                spectra_dict[date]['line'][line_name]['velocity'][param] = {}
                number_curves = len(spectra_dict[date]['line'][line_name]['line_extremes'][param])
                if None not in spectra_dict[date]['line'][line_name]['line_extremes'][param]:
                    wavelength_observed[param] = [spectra_dict[date]['line'][line_name]['line_extremes'][param][i]['x'] for i in range(number_curves)]
                    velocity_list = calculate_expansion_velocity(wavelength_expected, wavelength_observed[param])
                    spectra_dict[date]['line'][line_name]['velocity'][param]['mean'] = np.mean(velocity_list)
                    spectra_dict[date]['line'][line_name]['velocity'][param]['std'] = np.std(velocity_list)
                else:
                    spectra_dict[date]['line'][line_name]['velocity'][param]['mean'] = None
                    spectra_dict[date]['line'][line_name]['velocity'][param]['std'] = None
    return spectra_dict

def add_name_t_from_discovery_to_df(df, name, discovery_date):
    df['datetime'] = data_import.convert_to_mjd(df['datetime'], from_datetime=True)
    df['t_from_discovery'] = df['datetime'] - discovery_date
    df['name'] = name
    return df


def make_velocity_df(SN_dict, lines_dict):
    dates = SN_dict['spectra'].keys()
    lines = lines_dict.keys()
    velocity_df = []
    for date in dates:
        for line in lines:
            velocity_df.append(pd.DataFrame({'line': line,
                                             'absorption_mean_velocity': SN_dict['spectra'][date]['line'][line]['velocity']['absorption']['mean'],
                                             'absorption_std_velocity': SN_dict['spectra'][date]['line'][line]['velocity']['absorption']['std'],
                                             'emission_mean_velocity': SN_dict['spectra'][date]['line'][line]['velocity']['emission']['mean'],
                                             'emission_std_velocity':SN_dict['spectra'][date]['line'][line]['velocity']['emission']['std'],
                                             't_from_discovery': SN_dict['spectra'][date]['t_from_discovery']},
                                              index=[0]))
    velocity_df = pd.concat(velocity_df)
    velocity_df['name'] = SN_dict['name']
    return velocity_df


# plot expansion velocities over datetime
def plot_expansion_velocities(df_list, absorptionORemission):
    fig, ax = plt.subplots(1, figsize=(9, 6))

    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#d95f02'}
    markers_fill = {'SN2018hmx': 'full', 'SNiPTF14hls': 'none', 'SN2018aad': 'none'}
    names = []
    for df in df_list:
        SN_name = pd.unique(df['name'])[0]
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
    ax.set_title('Expansion velocity over time - '+absorptionORemission + ' - '+ str(names), fontsize=16)
    ax.set_ylabel('Expansion velocity (km/s)', size=12)
    ax.set_xlabel('Days from discovery', size=12)
    ax.legend()
    ax.tick_params(axis='both', which='major', labelsize=14)
    fig.savefig(''.join(names) + absorptionORemission + '_expansion_velocity_over_time_sn1999em' + '.png')

