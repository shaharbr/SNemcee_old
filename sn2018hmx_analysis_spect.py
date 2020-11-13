from matplotlib import pyplot as plt
import spectra_velocity
import data_import
import os
import pandas as pd
import re
import numpy as np

# plotting parameters
plt.rc('font', size=12)          # controls default text sizes
plt.rc('axes', titlesize=12)     # fontsize of the axes title
plt.rc('axes', labelsize=12)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=11)    # legend fontsize
plt.rc('figure', titlesize=14)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'

# For spectrum files taken from SNEx, get the mjd date of the observation from the filename
def mjd_from_filename(filename):
    date = re.sub('.*hmx|.*aad|_|-|P60.*|v1|[a-z]|[A-Z]|\..*', '', filename)
    # transform to standard datetime format
    date = pd.to_datetime(date)
    # convert date to MJD
    date_mjd = data_import.convert_to_mjd(date, from_datetime=True)
    return date_mjd

# For spectrum files taken from SNEx, get the telescope name from the filename
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

# add a spectrum to a dictionary that contains all the spectra of that SN. The spectrum is stored as a dataframe
# under the respective mjd date it was observed on.
def add_spectrum_to_dict(spectra_dict, filepath, telescope_name, date_mjd):
    spectra_dict[date_mjd] = {'telescope': telescope_name}
    # TODO add the dy of the spectrum to the uncertainties of the velocities, for the spectra recorded from ZTF?
    if telescope_name == 'P60':
        spectra_dict[date_mjd]['df'] = pd.read_csv(filepath, sep=' ', names=["x", "y", 'dy'], header=180)
    else:
        spectra_dict[date_mjd]['df'] = pd.read_csv(filepath, sep=' ', names=["x", "y"])
    return spectra_dict

# correct the time from discovery for redshift, by dividing the time by (1+z)
def add_rest_frame_days_from_discovery(SN_dict):
    dates = SN_dict['spectra'].keys()
    z = SN_dict['z']
    for date in dates:
        SN_dict['spectra'][date]['t_from_discovery'] = (date - SN_dict['discovery_date']) / (1 + z)
    return SN_dict

# correct wavelengths for redshift, by dividing the wavelengths by (1+z)
def correct_spectra_redshift(SN_dict):
    dates = SN_dict['spectra'].keys()
    z = SN_dict['z']
    for date in dates:
        SN_dict['spectra'][date]['df']['x'] = SN_dict['spectra'][date]['df']['x'] / (1 + z)
    return SN_dict

# normalize all the spectra by their mean, for display only
def normalize_spectra(SN_dict):
    dates = sorted(list(SN_dict['spectra'].keys()))
    for date in dates:
        avg = np.mean(SN_dict['spectra'][date]['df']['y'])
        SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / avg
        if avg > (10 ** -15):
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] * 2
        if avg > (10 ** -10):
            SN_dict['spectra'][date]['df']['y'] = SN_dict['spectra'][date]['df']['y'] / 4
    return SN_dict

# smooth the spectra proportionaly to their "noise" (std of the first 50 data points), for display only
def smooth_spectrum(spectrum_y, lineplot=False):
    var_list = list(spectrum_y)
    var_list = [x for x in var_list if str(x) != 'nan']
    var_list = var_list[-100:-1]
    start_variability = np.min([np.std(var_list) / 2, 0.3])
    if lineplot:
        smoothing_window = np.max([int(len(spectrum_y) / 70 * start_variability), 1])
    else:
        smoothing_window = np.max([int(len(spectrum_y) / 20 * start_variability), 1])
    smooth_y = spectrum_y.rolling(smoothing_window, center=True).mean()
    return smooth_y


def plot_stacked_spectra(SN_dict, lines_dict, plot_curve_fits=False, line_velocity=False):
    dates = sorted(SN_dict['spectra'].keys(), reverse=True)
    if line_velocity:
        fig, ax = plt.subplots(1, figsize=(3, 10))
    else:
        fig, ax = plt.subplots(1, figsize=(6, 10))
    name = SN_dict['Name']
    y_shift = 0.0
    # lines_list = lines_dict.keys()
    lines_list = ['Halpha', 'Hbeta', 'FeII 5169']
    colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#d95f02'}
    # draw the expected line wavelengths as vertical lines in background
    # if its a "lines velocity" type graph, draw the line at the x=0 (velocity=0) point
    if line_velocity:
        ax.axvline(x=0, color=colors[line_velocity], alpha=0.5)
    # if normal full spectra graph, draw the line in the respective wavelength of each line
    else:
        for line_name in lines_list:
            wavelength_expected = lines_dict[line_name]['peak']
            ax.axvline(x=wavelength_expected, color=colors[line_name], alpha=0.5)
    # draw the spectrum for each date, stacked with incremental shifts as defined by y_shift
    for date in dates:
        skip_date = False
        date_dict = SN_dict['spectra'][date]
        df = date_dict['df']
        # if its a "line velocity" type graph, draw only th e
        if line_velocity:
            # because of the blue-ness of the early spectra, any Pcygni other than Halpha
            # are very hard to identify (and therefore innacurate)
            if line_velocity != 'Halpha' and date_dict['t_from_discovery'] < 20:
                skip_date = True
            # smooth the spectrum for display, because its a line velocity graph it will be smoothed less
            smooth_y = smooth_spectrum(df['y'], line_velocity)
            # draw from 150 before to 200 after the expected Pcygni feature
            min_x = lines_dict[line_velocity]['absorption_range'][0] - 250
            max_x = lines_dict[line_velocity]['absorption_range'][1] + 250
            bool_x = (df['x'] > min_x) & (df['x'] < max_x)
            x = df['x'].loc[bool_x]
            # if the x is empty, its because the relevant segment of the spectra for this line doesn't exist
            if len(x) < 1:
                skip_date = True
            else:
                smooth_y = smooth_y[bool_x]
                # add the necessary shift for the stacking
                smooth_y = smooth_y + y_shift
                # replace the x values with the velocities, calculated based on the distance of the absorption trough
                # from the expected peak location of this line.
                wavelength_expected = lines_dict[line_velocity]['peak']
                x = spectra_velocity.calculate_expansion_velocity(wavelength_expected, x)
                # invert the x axis so the positive expansion velocities are to the left (like the shorter wavelegnths)
                ax.invert_xaxis()
                # plot title and labels
                ax.set_title(line_velocity)
                ax.set_xlabel('Velocity (km/s)')
                # make xlim tight around the min and max velocities
                ax.set_xlim(np.max(x), np.min(x))
                # define the location of the label for the telescope and rest-frame day of the observation
                label_x = np.min (x) - 150
                # define the shift for the next spectrum stacked on top of this one, based on the height of the present spectrum
                y_shift_delta = (np.max(smooth_y) - y_shift) / 10

        else:  # if its a normal stacked spectra plot
            smooth_y = smooth_spectrum(df['y'])
            smooth_y = smooth_y + y_shift
            x = df['x']
            # plot title, labels and xlim
            ax.set_title(name + ' spectra over time')
            ax.set_xlabel('Rest (z = 0.038) Wavelength (Å)')
            ax.set_xlim(3300, 9700)
            # define the location of the label for the telescope and rest-frame day of the observation
            label_x = 9780
            # define the shift for the next spectrum stacked on top of this one, based on the height of the present spectrum
            y_shift_delta = (np.max(smooth_y) - y_shift) / 2.5

        if not skip_date:
            # plot
            ax.plot(x, smooth_y, color='k', linewidth=1)
            # spectrum labeled with day from explosion
            spectrum_date_label = str(int(date_dict['t_from_discovery'])) + 'd (' + date_dict['telescope'] + ')'
            label_ypos = np.mean(smooth_y) - 0.7
            ax.text(label_x, label_ypos, spectrum_date_label, color='k', fontsize=10)
            # adding the polyfit curves if asked
            if plot_curve_fits:
                spectra_velocity.add_fitted_curves_to_plot(date_dict, lines_dict, ax, y_shift, line_velocity, number_curves_to_draw=100)
            y_shift += np.max([y_shift_delta, 0.3])
            if date_dict['t_from_discovery'] < 30:
                y_shift -= 0.5
            # if (date_dict['t_from_discovery'] > 80) & (date_dict['t_from_discovery'] < 83) :
            #     y_shift += 0.3
            # if date_dict['t_from_discovery'] >150:
            #     y_shift += 0.3
    ax.set_ylabel('Normalized fλ + shift')
    ax.set_yticks([])
    # ax.get_legend().remove()
    ax.tick_params(axis='both', which='major')
    fig.subplots_adjust(right=0.7)
    if plot_curve_fits:
        fig.savefig(r'figures/' + name + '_spectra_over_time_stacked_with_polyfit' + str(line_velocity) + '.png')
        fig.savefig(r'figures/' + name + '_spectra_over_time_stacked_with_polyfit' + str(line_velocity) + '.svg')
    else:
        fig.savefig(r'figures/' + name + '_spectra_over_time_stacked' + str(line_velocity) + '.png')
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
lines_dict = {'Halpha': {'peak': 6562.81, 'absorption_range': [6300, 6562], 'emission_range': [6500, 6700], 'width': 200},
              'Hbeta': {'peak': 4861, 'absorption_range': [4700, 4861], 'emission_range': [4800, 5000], 'width': 200},
              'FeII 5169': {'peak': 5169, 'absorption_range': [5080, 5169], 'emission_range': [5200, 5600], 'width': 150}}

# lines_dict = {'Halpha': {'peak': 6562.81, 'absorption_range': [6300, 6562], 'emission_range': [6500, 6700], 'width': 200},
#               'Hbeta': {'peak': 4861, 'absorption_range': [4700, 4861], 'emission_range': [4800, 5000], 'width': 200},
#               'FeII 5169': {'peak': 5169, 'absorption_range': [5080, 5169], 'emission_range': [5200, 5600], 'width': 150},
#               'FeII 5018': {'peak': 5018, 'absorption_range': [4930, 5010], 'emission_range': [5000, 5200], 'width': 70}}



# initialize sn18hmx dictionary
SN2018hmx = {'Name': 'SN2018hmx',
             'z': z_dict['SN2018hmx'],
             'discovery_date': discovery_date_dict['SN2018hmx'],
             'spectra': {}, 'expansion_velocities': []}

# add SNEx spectra to the 18hmx spectra dict
SNEX_dir = os.path.join('data', 'snexdata_target5025')
SNEX_filenames = list((file for file in os.listdir(SNEX_dir)
                      if os.path.isfile(os.path.join(SNEX_dir, file))))
for SNEXfile in SNEX_filenames:
    SNEX_filepath = os.path.join(SNEX_dir, SNEXfile)
    SN2018hmx['spectra'] = add_spectrum_to_dict(SN2018hmx['spectra'],
                                     SNEX_filepath,
                                     telescope_from_filename(SNEX_filepath),
                                     mjd_from_filename(SNEX_filepath))
# add Keck spectra to the 18hmx spectra dict
keck_spectra_dir = os.path.join('data', 'uncalibrated_keck_spectra')
SN2018hmx['spectra'] = add_spectrum_to_dict(SN2018hmx['spectra'],
                                 os.path.join(keck_spectra_dir, '2018hmx_20190111_2458494.90432_1.ascii'),
                                 'Keck',
                                 data_import.convert_to_mjd(2458494.90432))
SN2018hmx['spectra'] = add_spectrum_to_dict(SN2018hmx['spectra'],
                                 os.path.join(keck_spectra_dir, '2018hmx_20191023_2458780.09736_1.ascii'),
                                 'Keck',
                                 data_import.convert_to_mjd(2458780.09736))

# save_spect_csv(SN18hmx_spect_dict, 'snexdata_target5025')


# initialize sn18aad dictionary
SN2018aad = {'Name': 'SN2018aad',
             'z': z_dict['SN2018aad'],
             'discovery_date': discovery_date_dict['SN2018aad'],
             'spectra': {}, 'expansion_velocities': []}

# add SNEx spectra to the 18aad spectra dict
SNEX_dir = os.path.join('data', 'snexdata_target4771')
SNEX_filenames = list((file for file in os.listdir(SNEX_dir)
                      if os.path.isfile(os.path.join(SNEX_dir, file))))
for SNEXfile in SNEX_filenames:
    SNEX_filepath = os.path.join(SNEX_dir, SNEXfile)
    SN2018aad['spectra'] = add_spectrum_to_dict(SN2018aad['spectra'],
                                     SNEX_filepath,
                                     telescope_from_filename(SNEX_filepath),
                                     mjd_from_filename(SNEX_filepath))

# initialize sn18aad dictionary
SNiPTF14hls = {'Name': 'SNiPTF14hls',
               'z': z_dict['SNiPTF14hls'],
               'discovery_date': discovery_date_dict['SNiPTF14hls'],
               'spectra': {},
               'expansion_velocities': []}

# import existing iPTF14hls expansion velocity file to the dict of iPTF14hls
SN14hls_veloc_path = os.path.join('results', 'iPTF14hls_expansion_velocity.csv')
SN14hls_veloc_df = pd.read_csv(SN14hls_veloc_path, header=0)
SN14hls_veloc_df['JD'] = pd.to_datetime(SN14hls_veloc_df['JD'], unit='D', origin='julian')
SN14hls_veloc_df.rename(columns={'JD': 'datetime', 'Velocity [km/s]': 'absorption_mean_velocity', 'Line': 'line',
                                       'Velocity_Error [km/s]': 'absorption_std_velocity'}, inplace=True)
SNiPTF14hls['expansion_velocities'] = SN14hls_veloc_df



# TODO fix bugs in pEW

# for each SNe: add restframe time, redshift corrections, normalization and plotting
for SN in [SN2018hmx, SN2018aad]:
    # add rest frame days from discovery in 't_from_disovery'
    SN = add_rest_frame_days_from_discovery(SN)
    # correct spectra for redshift
    SN = correct_spectra_redshift(SN)
    # normalize spectra for visualization
    SN = normalize_spectra(SN)
    # plot stacked spectra
    plot_stacked_spectra(SN, lines_dict)
    # produce Pcygni curve fits and add them to the SN dict
    SN['spectra'] = spectra_velocity.fit_Pcygni_curves(SN['spectra'], lines_dict, fixed_curve_range=False, number_curves=100)
    # plot stacked spectra with the Pcygni curves found
    plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True)
    # make a reduced dictionary with only the first 170 days (a bit more than when we can expect real Pcygni features)

    print(len(SN))
    print(SN['spectra'].keys())
    SN_170d = {'Name': SN['Name'], 'spectra': {}, 'expansion_velocities': []}
    dates = sorted(list(SN['spectra'].keys()))
    for date in dates:
        if int(SN['spectra'][date]['t_from_discovery']) < 170:
            print(SN['spectra'][date]['t_from_discovery'])
            SN_170d['spectra'][date] = SN['spectra'][date]
    print(len(SN_170d))
    print(SN_170d['spectra'].keys())
    # produce plots showing in zoom-in the Pcygni for each line, demonstrating the velocities
    plot_stacked_spectra(SN_170d, lines_dict, plot_curve_fits=True, line_velocity='Halpha')
    plot_stacked_spectra(SN_170d, lines_dict, plot_curve_fits=True, line_velocity='Hbeta')
    plot_stacked_spectra(SN_170d, lines_dict, plot_curve_fits=True, line_velocity='FeII 5169')
    # calculate expansion velocities from the Pcygni fits found
    SN_170d['spectra'] = spectra_velocity.add_expansion_velocity(SN_170d['spectra'], lines_dict)
    # make a dataframe summerizing the velocities
    SN_170d['expansion_velocities'] = spectra_velocity.make_velocity_df(SN_170d, lines_dict)
    # plot_pEW(SN)
    spectra_velocity.plot_expansion_velocities([SN_170d['expansion_velocities']], 'absorption')
    SN_170d['expansion_velocities'].to_csv(os.path.join('results', SN_170d['Name']+'_expansion_velocities.csv'))

# plot expansion velocities of SNe against each other
# print(SN2018hmx['expansion_velocities'])
# print(SNiPTF14hls['expansion_velocities'])
# spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities']], 'absorption')
# TODO fix the way 14hls velocities are imported, something isn't working there
# spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SNiPTF14hls['expansion_velocities']], 'absorption')
# spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SN2018aad['expansion_velocities']], 'absorption')
# spectra_velocity.plot_expansion_velocities([SN2018aad['expansion_velocities']], 'absorption')

# save expansion velocities to csv
# SN2018hmx['expansion_velocities'].to_csv(os.path.join('results', 'sN2018hmx_expansion_velocities.csv'))
# SNiPTF14hls['expansion_velocities'].to_csv(os.path.join('results', 'SNiPTF14hls_expansion_velocities.csv'))
# SN2018aad['expansion_velocities'].to_csv(os.path.join('results', 'SN2018aad_expansion_velocities.csv'))

plt.show()

