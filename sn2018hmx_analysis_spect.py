from matplotlib import pyplot as plt
import spectra_velocity
import data_import
import os
import pandas as pd
import re
import numpy as np
# plotting parameters
plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'


def mjd_from_filename(filename):
    date = re.sub('.*hmx|.*aad|_|-|P60.*|v1|[a-z]|[A-Z]|\..*', '', filename)
    # transform to standard datetime format
    date = pd.to_datetime(date)
    # convert date to MJD
    date_mjd = data_import.convert_to_mjd(date, from_datetime=True)
    return date_mjd

def telescope_from_filename(filename):
    if 'ZTF' in os.path.basename(filename):
        telescope_name = 'P60'
    elif 'redblu' in os.path.basename(filename):
        telescope_name = 'Las Cumbres'
    elif 'HET' in os.path.basename(filename):
        telescope_name = 'HET'
    else:
        telescope_name = 'ND'
    return telescope_name

def add_spectrum_to_dict(spectra_dict, filepath, telescope_name, date_mjd):
    spectra_dict[date_mjd] = {'telescope': telescope_name}
    if telescope_name == 'ZTF':
        spectra_dict[date_mjd]['df'] = pd.read_csv(filepath, sep=' ', names=["x", "y", 'dy'], header=180)
    else:
        spectra_dict[date_mjd]['df'] = pd.read_csv(filepath, sep=' ', names=["x", "y"])
    return spectra_dict

'''
# TODO adapt this func
def make_SN_dict(SN_name, z, discovery_date_mjd, spectra_dict):
    # organize in SN_dict
    fields = [lightcurves_dict, z_dict, discovery_date_dict,
                 distance_modulus_dict, galactic_extinction_dict,
                 spectra, expansion_velocities]
    keys = ['lightcurve', 'z','discovery_date', 'distance_modulus', 'galactic_extinction',
               'spectra', 'expansion_velocities']

    # convert date to MJD
    discovery_date_dict[SN_name] = convert_to_mjd(discovery_date_dict[SN_name], from_datetime=True)
    SN_dict = {}
    SN_dict['Name'] = SN_name
    for i in range(len(fields)):
        if fields[i]:
            SN_dict[keys[i]] = fields[i][SN_name]
        else:
            SN_dict[keys[i]] = ''
    return SN_dict
'''


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
        avg = np.mean(SN_dict['spectra'][date]['df']['y'])
        SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / avg
        if SN_dict['spectra'][date]['telescope'] == 'P60':
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] * 2
        if SN_dict['spectra'][date]['telescope'] == 'HET':
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / 2
    return SN_dict


def smooth_spectrum(spectrum_y, lineplot=False):
    var_list = list(spectrum_y)
    var_list = [x for x in var_list if str(x) != 'nan']
    var_list = var_list[0:100]
    start_variability = np.min([np.std(var_list) / 2, 0.2])
    if lineplot:
        smoothing_window = int(len(spectrum_y) / 70 * start_variability) + 1
    else:
        smoothing_window = int(len(spectrum_y)/30 * start_variability)+1
    smooth_y = spectrum_y.rolling(smoothing_window, center=True).mean()
    return smooth_y


# plot all spectra in an overlay figure
def plot_overlay_spectra(SN_dict, lines_dict=False):
    dates = SN_dict['spectra'].keys()
    fig, ax = plt.subplots(1, figsize=(10, 5))
    name = SN_dict['Name']
    for date in dates:
        SN_dict['spectra'][date]['df'].plot(x='x', y='y', ax=ax, label=date)
    if lines_dict:
        for line_name in lines_dict.keys():
            ax.axvline(x=lines_dict[line_name]['peak'], color='k')
    ax.set_title(name+' spectra over time')
    ax.set_xlabel('Rest (z = 0.038) Wavelength (Å)')
    ax.set_ylabel('Flux (10-15 erg s-1 cm-2 Å-1)')
    ax.tick_params(axis='both', which='major')
    fig.savefig(r'figures/' +name + '_spectra_over_time_overlay'+'.png')
    fig.savefig(r'figures/' + name + '_spectra_over_time_overlay' + '.svg')


def plot_stacked_spectra(SN_dict, lines_dict, plot_curve_fits=False, line_velocity=False):
    dates = sorted(SN_dict['spectra'].keys(), reverse=True)
    if line_velocity:
        fig, ax = plt.subplots(1, figsize=(5, 10))
    else:
        fig, ax = plt.subplots(1, figsize=(8, 8))
    name = SN_dict['Name']
    y_shift = 0.0
    y_shift_delta = 0.0
    lines_list = lines_dict.keys()
    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#d95f02', 'FeII 5018': '#1D78EF'}
    # draw expected lines locations in background
    if line_velocity:
        ax.axvline(x=0, color=colors[line_velocity], alpha=0.5)
    else:
        for line_name in lines_list:
            wavelength_expected = lines_dict[line_name]['peak']
            # TODO remove 5018 line
            ax.axvline(x=wavelength_expected, color=colors[line_name], alpha=0.5)
        ax.axvline(x=5018, color='#1D78EF', alpha=0.5)
    # draw spectra
    for date in dates:
        date_dict = SN_dict['spectra'][date]
        df = date_dict['df']
        # if date_dict['t_from_discovery'] > 300:
        #     print('late')
        #     y_shift_delta = np.max([y_shift_delta, 2])

        # if date_dict['t_from_discovery'] < 20:
        #     print()
        # else:
        # y_shift += y_shift_delta
        y_shift += np.max([y_shift_delta, 1])
        if line_velocity:
            smooth_y = smooth_spectrum(df['y'], lineplot=True)
            # if line_velocity != 'Halpha':
            #     y_shift_delta = 0.3
            #     if date_dict['t_from_discovery'] < 20:
            #         y_shift = y_shift - 0.6
            # else:
            #     y_shift_delta = 0.8
            # y_shift_delta = 1
            min_x = lines_dict[line_velocity]['absorption_range'][0] - 150
            if date_dict['t_from_discovery'] < 70:
                min_x = min_x - 100
            max_x = lines_dict[line_velocity]['peak'] + 200
            bool_x = (df['x'] > min_x) & (df['x'] < max_x)
            x = df['x'].loc[bool_x]
            y = smooth_y[bool_x]
            y_shift_delta = (np.max(y) - np.min(y)) / 3
            # y = df['y'].loc[bool_x]
            wavelength_expected = lines_dict[line_velocity]['peak']
            x = spectra_velocity.calculate_expansion_velocity(wavelength_expected, x)
            if len(x) > 0:
                max_x = np.nanmax(x)
                min_x = np.nanmin(x)
                ax.set_xlim(min_x, max_x)
                label_x = min_x - 150
            y = y + y_shift
            ax.invert_xaxis()
            ax.set_title(line_velocity)
            ax.set_xlabel('Velocity (km/s)')
        else:
            smooth_y = smooth_spectrum(df['y'])
            x = df['x']
            # y_shift_delta = (np.max(smooth_y) - np.min(smooth_y)) / 3
            y_shift_delta = (np.max(smooth_y) - np.min(smooth_y)) / 2.5
            y = smooth_y + y_shift
            ax.set_xlim(3300, 9700)
            label_x = 9780
            ax.set_title(name + ' spectra over time')
            ax.set_xlabel('Rest (z = 0.038) Wavelength (Å)')
        if len(x) > 0 and len(y) > 0:
            # plot
            ax.plot(x, y, color='k', linewidth=1)
            # spectrum labeled with day from explosion
            spectrum_date_label = str(int(date_dict['t_from_discovery'])) + 'd (' + date_dict['telescope'] + ')'
            label_ypos = np.mean(y) - 0.7
            ax.text(label_x, label_ypos, spectrum_date_label, color='k', fontsize=10)
            # adding the polyfit curves if asked
            if plot_curve_fits:
                spectra_velocity.add_fitted_curves_to_plot(date_dict, lines_dict, ax, y_shift, line_velocity, number_curves_to_draw=10)
        # y_shift += y_shift_delta
    ax.set_ylabel('Normalized fλ + shift')
    ax.set_yticks([])
    # ax.get_legend().remove()
    ax.tick_params(axis='both', which='major')
    fig.subplots_adjust(right=0.7)
    if plot_curve_fits:
        fig.savefig(r'figures/' + name + '_spectra_over_time_stacked_with_polyfit' + str(line_velocity) + '.png')
        fig.savefig(r'figures/' + name + '_spectra_over_time_stacked_with_polyfit' + str(line_velocity) + '.svg')

    else:
        fig.savefig(r'figures/' +name + '_spectra_over_time_stacked' + str(line_velocity) + '.png')
        fig.savefig(r'figures/' + name + '_spectra_over_time_stacked' + str(line_velocity) + '.svg')


def plot_pEW(SN):
    dates = SN['spectra'].keys()
    pEW_SN = [[], [], []]
    for date in sorted(dates):
        result, result_std = spectra_velocity.pEW(SN['spectra'], 'FeII 5018', date)
        pEW_SN[0].append(SN['spectra'][date]['t_from_discovery'])
        pEW_SN[1].append(result)
        pEW_SN[2].append(result_std)
    plt.figure()
    plt.title('FeII 5018 pEW over time')
    pEW_plot = [[], [], []]
    pEW_plot[0] = [pEW_SN[0][i] for i in [3, 4, 5, 6, 7, 9, 10, 11]]
    pEW_plot[1] = [pEW_SN[1][i] for i in [3, 4, 5, 6, 7, 9, 10, 11]]
    pEW_plot[2] = [pEW_SN[2][i] for i in [3, 4, 5, 6, 7, 9, 10, 11]]
    plt.plot(pEW_plot[0], pEW_plot[1])
    plt.errorbar(x=pEW_plot[0], y=pEW_plot[1], yerr=pEW_plot[2], marker='o', linestyle='None', )
    plt.ylim(0, 40)
    plt.xlim(0, 120)
    plt.xlabel('days from discovery')
    plt.ylabel('pEW')

def save_spect_csv(dict, dir):
    names = dict.keys()
    os.mkdir(os.path.join('data', dir, 'csv'))
    # for name in names:
        # dict[name]['df'].to_csv(os.path.join('data', dir, 'csv', str(name)+'.csv'))


# dictionary of redshifts for each SNe
z_dict = {'SN2018hmx': 0.038,
     'SN1999em': 0.0024,
     'SNiPTF14hls': 0.0344,
     'SN2018aad': 0.023}

# dictionary of discovery dates for each SNe
# TODO make sure 14hls time from discovery is true
discovery_date_dict = {'SN2018hmx': 58408.64599537,
                  'SN1999em': 51480.43958333,
                  'SNiPTF14hls': 56922.52986111,
                  'SN2018aad': 58182.79959491}

# dictionary describing for each line (Halpha, Hbeta, FeII 5169 and FeII 5018), where its expected peak should be,
# the range where the emission should be found, the range where the blueshifted absotbtion troughs are expected to
# be found, and the expected width of the curves that should be fitted to the peaks and throughs of the
# P Cygni features. Determined uniformly for all spectra and SNe, based on visual inspection of the spectra.
lines_dict = {'Halpha': {'peak': 6562.81, 'absorption_range': [6200, 6562], 'emission_range': [6500, 6700], 'width': 200},
              'Hbeta': {'peak': 4861, 'absorption_range': [4700, 4861], 'emission_range': [4800, 5000], 'width': 200},
              'FeII 5169': {'peak': 5169, 'absorption_range': [5080, 5169], 'emission_range': [5200, 5600], 'width': 150},
              'FeII 5018': {'peak': 5018, 'absorption_range': [4930, 5010], 'emission_range': [5000, 5200], 'width': 70}}

# initialize sn18hmx dictionary
SN2018hmx = {'z': z_dict['SN2018hmx'],
             'discovery_date': discovery_date_dict['SN2018hmx'],
             'spectra': {}, 'expansion_velocities': []}

# add SNEx spectra to the 18hmx spectra dict
SNEX_dir = os.path.join('data', 'snexdata_target5025')
SNEX_filenames = list((file for file in os.listdir(SNEX_dir)
                      if os.path.isfile(os.path.join(SNEX_dir, file))))
for SNEXfile in SNEX_filenames:
    SNEX_filepath = os.path.join(SNEX_dir, SNEXfile)
    SN2018hmx = add_spectrum_to_dict(SN2018hmx,
                                     SNEX_filepath,
                                     telescope_from_filename(SNEX_filepath),
                                     mjd_from_filename(SNEX_filepath))
# add Keck spectra to the 18hmx spectra dict
keck_spectra_dir = os.path.join('data', 'uncalibrated_keck_spectra')
SN2018hmx = add_spectrum_to_dict(SN2018hmx,
                                 os.path.join(keck_spectra_dir, '2018hmx_20190111_2458494.90432_1.ascii'),
                                 'Keck',
                                 data_import.convert_to_mjd(2458494.90432))
SN2018hmx = add_spectrum_to_dict(SN2018hmx,
                                 os.path.join(keck_spectra_dir, '2018hmx_20191023_2458780.09736_1.ascii'),
                                 'Keck',
                                 data_import.convert_to_mjd(2458780.09736))

# save_spect_csv(SN18hmx_spect_dict, 'snexdata_target5025')


# initialize sn18aad dictionary
SN2018aad = {'z': z_dict['SN2018aad'],
             'discovery_date': discovery_date_dict['SN2018aad'],
             'spectra': {}, 'expansion_velocities': []}

# add SNEx spectra to the 18aad spectra dict
SNEX_dir = os.path.join('data', 'snexdata_target4771')
SNEX_filenames = list((file for file in os.listdir(SNEX_dir)
                      if os.path.isfile(os.path.join(SNEX_dir, file))))
for SNEXfile in SNEX_filenames:
    SNEX_filepath = os.path.join(SNEX_dir, SNEXfile)
    SN2018aad = add_spectrum_to_dict(SN2018aad,
                                     SNEX_filepath,
                                     telescope_from_filename(SNEX_filepath),
                                     mjd_from_filename(SNEX_filepath))


# spectra = {'SN2018hmx': SN2018aad}
           # 'SN2018aad': sn18aad_spectra_dict}


# SN2018hmx = data_import.make_SN_dict('SN2018hmx', z_dict=z, discovery_date_dict=discovery_date, spectra=spectra)
# SN2018aad = data_import.make_SN_dict('SN2018aad', z_dict=z, discovery_date_dict=discovery_date, spectra=spectra)
# SNiPTF14hls = data_import.make_SN_dict('SNiPTF14hls', z_dict=z, discovery_date_dict=discovery_date, expansion_velocities=expansion_v)


# import existing expansion velocity file for iPTF14hls
SN14hls_veloc_path = os.path.join('results', 'iPTF14hls_expansion_velocity.csv')
SN14hls_veloc_df = pd.read_csv(SN14hls_veloc_path, header=0)
SN14hls_veloc_df['JD'] = pd.to_datetime(SN14hls_veloc_df['JD'], unit='D', origin='julian')
SN14hls_veloc_df.rename(columns={'JD': 'datetime', 'Velocity [km/s]': 'absorption_mean_velocity', 'Line': 'line',
                                       'Velocity_Error [km/s]': 'absorption_std_velocity'}, inplace=True)

# initialize sn18aad dictionary
SNiPTF14hls = {'z': z_dict['SNiPTF14hls'],
               'discovery_date': discovery_date_dict['SNiPTF14hls'],
               'spectra': {},
               'expansion_velocities': SN14hls_veloc_df}


# expansion_v = {'SNiPTF14hls': SN14hls_expans_v_df}
# SNiPTF14hls['expansion_velocities'] = spectra_velocity.add_name_t_from_discovery_to_df(SNiPTF14hls['expansion_velocities'], 'SNiPTF14hls', discovery_date['SNiPTF14hls'])

# TODO fix bugs in pEW

# for each SNe: add restframe time, redshift corrections, normalization and plotting
for SN in [SN2018hmx]:
    # add rest frame days from discovery in 't_from_disovery'
    SN = add_rest_frame_days_from_discovery(SN)
    # correct spectra for redshift
    SN = correct_specta_redshift(SN)
    # normalize spectra for visualization
    SN = normalize_spectra(SN)
    # plot stacked spectra
    plot_stacked_spectra(SN, lines_dict)
    # produce Pcygni curve fits and add them to the SN dict
    SN['spectra'] = spectra_velocity.fit_Pcygni_curves(SN['spectra'], lines_dict, fixed_curve_range=False, number_curves=60)
    # plot stacked spectra with the Pcygni curves found
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True)
    # produce plots showing in zoom-in the Pcygni for each line, demonstrating the velocities
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='Halpha')
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='Hbeta')
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='FeII 5169')
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='FeII 5018')
    # calculate expansion velocities from the Pcygni fits found
    SN['spectra'] = spectra_velocity.add_expansion_velocity(SN['spectra'], lines_dict)
    # make a dataframe summerizing the velocities
    SN['expansion_velocities'] = spectra_velocity.make_velocity_df(SN, lines_dict)
    # plot_pEW(SN)

# plot expansion velocities of SNe against each other
spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SNiPTF14hls['expansion_velocities']], 'absorption')
# spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SN2018aad['expansion_velocities']], 'absorption')

# save expansion velocities to csv
SN2018hmx['expansion_velocities'].to_csv(os.path.join('results', 'sN2018hmx_expansion_velocities.csv'))
SNiPTF14hls['expansion_velocities'].to_csv(os.path.join('results', 'SNiPTF14hls_expansion_velocities.csv'))
# SN2018aad['expansion_velocities'].to_csv(os.path.join('results', 'SN2018aad_expansion_velocities.csv'))

plt.show()

