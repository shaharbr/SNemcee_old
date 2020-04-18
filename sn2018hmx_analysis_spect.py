from matplotlib import pyplot as plt
import pandas as pd
import spectra_velocity
import data_import

#  TODO add the import from csv for all correction parameters like in the photometry analysis code

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



z = {'SN2018hmx': 0.038,
     'SN1999em': 0.0024,
     'SNiPTF14hls': 0.0344,
     'SN2018aad': 0.023}

discovery_date = {'SN2018hmx': '2018-10-17 15:30:14',
                  'SN1999em': '1999-10-29 10:33:00',
                  'SNiPTF14hls': '2014-09-22 12:43:00',
                  'SN2018aad': '2018-03-05 19:11:25'}

lines_dict = {'Halpha': {'peak': 6562.81, 'absorption_range': [6150, 6500], 'emission_range': [6400, 6600], 'width': 180},
              'Hbeta': {'peak': 4861, 'absorption_range': [4600, 4800], 'emission_range': [4700, 4900], 'width': 100},
              'FeII 5169': {'peak': 5169, 'absorption_range': [4950, 5150], 'emission_range': [5100, 5500], 'width': 120},
              'FeII 5018': {'peak': 5018, 'absorption_range': [4870, 4970], 'emission_range': [4900, 5100], 'width': 60}}


SN18hmx_spect_dict = data_import.LCO_spect(r'data\snexdata_target5025')
SN18aad_spect_dict = data_import.LCO_spect(r'data\snexdata_target4771')

# import expansion velocity file for iPTF14hls
SN14hls_expans_v_df = data_import.SN14hls_expans_v(r'results\iPTF14hls_expansion_velocity.csv')
# TODO correct 14hls time from discovery

spectra = {'SN2018hmx': SN18hmx_spect_dict,
           'SN2018aad': SN18aad_spect_dict}

expansion_v = {'SNiPTF14hls': SN14hls_expans_v_df}


SN2018hmx = data_import.make_SN_dict('SN2018hmx', z_dict=z, discovery_date_dict=discovery_date, spectra=spectra)
SN2018aad = data_import.make_SN_dict('SN2018aad', z_dict=z, discovery_date_dict=discovery_date, spectra=spectra)
SNiPTF14hls = data_import.make_SN_dict('SNiPTF14hls', z_dict=z, discovery_date_dict=discovery_date, expansion_velocities=expansion_v)

SNiPTF14hls['expansion_velocities'] = spectra_velocity.add_name_t_from_discovery_to_df(SNiPTF14hls['expansion_velocities'], 'SNiPTF14hls', discovery_date['SNiPTF14hls'])

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



for SN in [SN2018hmx, SN2018aad]:
    SN = spectra_velocity.add_rest_frame_days_from_discovery(SN)
    SN = spectra_velocity.correct_specta_redshift(SN)
    SN = spectra_velocity.normalize_spectra(SN)
    # TODO fix bugs in pEW

    spectra_velocity.plot_stacked_spectra(SN, lines_dict)
    SN['spectra'] = spectra_velocity.fit_Pcygni_curves(SN['spectra'], lines_dict, fixed_curve_range=False, number_curves=30)
    # spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True)
    # spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='Halpha')
    # spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='Hbeta')
    # spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='FeII 5169')
    # spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='FeII 5018')
    # SN['spectra'] = spectra_velocity.add_expansion_velocity(SN['spectra'], lines_dict)
    # SN['expansion_velocities'] = spectra_velocity.make_velocity_df(SN, lines_dict)
    plot_pEW(SN)


# spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SNiPTF14hls['expansion_velocities']], 'absorption')
# spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SN2018aad['expansion_velocities']], 'absorption')

# SN2018hmx['expansion_velocities'].to_csv(r'results\sN2018hmx_expansion_velocities.csv')
# SNiPTF14hls['expansion_velocities'].to_csv(r'results\SNiPTF14hls_expansion_velocities.csv')
# SN2018aad['expansion_velocities'].to_csv(r'results\SN2018aad_expansion_velocities.csv')

plt.show()

