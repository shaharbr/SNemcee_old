import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
import re
import random

# TODO: convert all datetime to MJD

plt.rcParams['font.sans-serif'] = 'Arial'



def add_time_from_discovery(SN_dict):
    dates = SN_dict['spectra'].keys()
    for date in dates:
        SN_dict['spectra'][date]['t_from_discovery'] = date - SN_dict['discovery_date']
    return SN_dict


def correct_redshift(SN_dict):
    dates = SN_dict['spectra'].keys()
    for date in dates:
        # correct wavelengths for redshift
        SN_dict['spectra'][date]['df']['x'] = SN_dict['spectra'][date]['df']['x'] / (1 + SN_dict['z'])
    return SN_dict

def normalize_spectra(SN_dict):
    dates = SN_dict['spectra'].keys()
    for date in dates:
        avg = np.average(SN_dict['spectra'][date]['df']['y'])
        SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / avg
        if SN_dict['spectra'][date]['telescope'] == 'HET':
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / 2
    return SN_dict

def smooth_LCO_spectra(SN_dict):
    dates = SN_dict['spectra'].keys()
    for date in dates:
        if SN_dict['spectra'][date]['telescope'] == 'LCO':
            start_variability = np.std(SN_dict['spectra'][date]['df']['y'][0:300])
            smoothing_window = int(50 * start_variability)
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'].rolling(smoothing_window).mean()
    return SN_dict


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
            slice_poly = np.polyfit(df_slice['x'], df_slice['y'], deg=2)
            f[param][i] = (np.poly1d(slice_poly))
            x_slice[param][i] = (df_slice['x'])
            curve_extreme[param][i] = find_curve_extreme(f[param][i], x_slice[param][i], param)
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


def plot_stacked_spectra(SN_dict, lines_dict=False, plot_curve_fits=False):
    dates = sorted(SN_dict['spectra'].keys(), reverse=True)
    fig, ax = plt.subplots(1, figsize=(7, 9))
    name = SN_dict['name']
    y_shift = 0
    y_shift_delta = 1.5
    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#''d95f02'}
    smoothed_SN_dict = smooth_LCO_spectra(SN_dict)
    if lines_dict:
        for line_name in lines_dict.keys():
            ax.axvline(x=lines_dict[line_name]['peak'], color=colors[line_name], alpha=0.5)
    for date in dates:
        smoothed_SN_dict['spectra'][date]['df']['y'] += y_shift
        smoothed_SN_dict['spectra'][date]['df'].plot(x='x', y='y', ax=ax, color='navy', linewidth=1)
        timedelta = re.sub('days.*', 'd', str(smoothed_SN_dict['spectra'][date]['t_from_discovery']))
        spectrum_label = timedelta + ' (' + smoothed_SN_dict['spectra'][date]['telescope'] + ')'
        label_ypos = np.mean(smoothed_SN_dict['spectra'][date]['df']['y']) - 0.4
        ax.text(9780, label_ypos, spectrum_label , fontsize=14, color='navy', fontname="Arial")
        if plot_curve_fits:
            add_fitted_curves_to_plot(smoothed_SN_dict['spectra'][date], lines_dict, ax, y_shift, number_curves_to_draw=3)
        y_shift += y_shift_delta
    ax.set_title(name + ' spectra over time', fontsize=16)
    ax.set_xlabel('Rest (z = 0.037) Wavelength (Å)', size=14)
    ax.set_ylabel('Normalized fλ + constant', size=14)
    ax.set_xlim(3300,9700)
    ax.set_yticks([])
    ax.get_legend().remove()
    ax.tick_params(axis='both', which='major', labelsize=14)
    fig.subplots_adjust(right=0.8)
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles[::-1], labels[::-1], title='Time after discovery', loc='upper right')
    # if plot_curve_fits:
        # fig.savefig(name + '_spectra_over_time_stacked_with_polyfit.png')
    # else:
        # fig.savefig(name + '_spectra_over_time_stacked.png')



def add_fitted_curves_to_plot(single_spectrum_dict, lines_dict, ax, y_shift=False, number_curves_to_draw=1):
    colors = ['r', 'g', 'b']
    n = 0
    for line_name in lines_dict.keys():
        color = colors[n]
        n +=1
        for param in ['absorption', 'emission']:
            for i in range(number_curves_to_draw):
                x_slice = single_spectrum_dict['line'][line_name]['line_xslices'][param][i]
                f = single_spectrum_dict['line'][line_name]['polyfits'][param][i]
                curve_y = f(x_slice)
                extreme_point = single_spectrum_dict['line'][line_name]['line_extremes'][param][i]
                if y_shift:
                    # TODO understand and fix this bug - why do I need to add y_shift twice here?
                    curve_y += y_shift
                    curve_y += y_shift
                    extreme_point['y'] += y_shift
                    extreme_point['y'] += y_shift
                ax.plot(x_slice, curve_y, color=color, alpha=0.5, linewidth=1.5)
                ax.plot(extreme_point['x'], extreme_point['y'], '.', color=color)


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
                wavelength_observed[param] = [spectra_dict[date]['line'][line_name]['line_extremes'][param][i]['x'] for i in range(number_curves)]
                velocity_list = calculate_expansion_velocity(wavelength_expected, wavelength_observed[param])
                spectra_dict[date]['line'][line_name]['velocity'][param]['mean'] = np.mean(velocity_list)
                spectra_dict[date]['line'][line_name]['velocity'][param]['std'] = np.std(velocity_list)
    return spectra_dict


def make_velocity_df(SN_dict, lines_dict):
    dates = SN_dict['spectra'].keys()
    lines = lines_dict.keys()
    velocity_df = []
    for date in dates:
        for line in lines:
            velocity_df.append(pd.DataFrame({'datetime': date,
                                             'line': line,
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

    names = []
    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#d95f02'}
    markers_fill = {'2018hmx': 'full', 'iPTF14hls': 'none'}
    for df in df_list:
        lines = colors.keys()
        SN_name = df['name'].unique()[0]
        names.append(SN_name)
        for line_name in lines:
            linedf = df.loc[df['line'] == line_name]
            timedelta = linedf.loc[linedf['line'] == line_name]['t_from_discovery'].dt.days

            ax.errorbar(x=timedelta, y=linedf[absorptionORemission+'_mean_velocity'],
                          yerr=linedf[absorptionORemission+'_std_velocity'],
                                                 label=SN_name+' '+line_name,
                                                 marker='s',
                                                 fillstyle=markers_fill[SN_name],
                                                 linestyle='None',
                                                 color=colors[line_name])




    if absorptionORemission == 'absorption':
        ax.set_ylim(bottom=0)
    ax.set_title('Expansion velocity over time - '+absorptionORemission, fontsize=16)
    ax.set_ylabel('Expansion velocity (km/s)', size=12)
    ax.set_xlabel('Days from discovery', size=12)
    ax.legend()
    ax.tick_params(axis='both', which='major', labelsize=14)
    fig.savefig(''.join(names) + absorptionORemission + '_expansion_velocity_over_time_sn1999em' + '.png')

