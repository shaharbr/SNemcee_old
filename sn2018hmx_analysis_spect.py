from matplotlib import pyplot as plt
import pandas as pd
import spectra_velocity
import data_import

plt.rcParams['font.sans-serif'] = 'Arial'


plt.rcParams['font.sans-serif'] = 'Arial'


z = {'SN2018hmx': 0.037,
     'SN1999em': 0.0024,
     'SNiPTF14hls': 0.0344,
     'SN2018aad': 0.025}

discovery_date = {'SN2018hmx': '2018-10-17 15:30:14',
                  'SN1999em': '1999-10-29 10:33:00',
                  'SNiPTF14hls': '2014-09-22 12:43:00',
                  'SN2018aad': '2018-03-05 01:55:12'} # TODO find real discovery date of 18aad

lines_dict = {'Halpha': {'peak': 6562.81, 'absorption_range': [6200, 6500], 'emission_range': [6400, 6600]},
              'Hbeta': {'peak': 4861, 'absorption_range': [4600, 4800], 'emission_range': [4700, 4900]},
              'FeII 5169': {'peak': 5169, 'absorption_range': [5000, 5100], 'emission_range': [5100, 5500]}}


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

# for SN in [SN2018hmx, SN2018aad]:
for SN in [SN2018hmx]:
    SN = spectra_velocity.add_rest_frame_days_from_discovery(SN)
    SN = spectra_velocity.correct_specta_redshift(SN)
    SN = spectra_velocity.normalize_spectra(SN)
    spectra_velocity.plot_stacked_spectra(SN, lines_dict)
    SN['spectra'] = spectra_velocity.fit_Pcygni_curves(SN['spectra'], lines_dict, fixed_curve_range=False, number_curves=10)
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True)
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='Halpha')
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='Hbeta')
    spectra_velocity.plot_stacked_spectra(SN, lines_dict, plot_curve_fits=True, line_velocity='FeII 5169')
    SN['spectra'] = spectra_velocity.add_expansion_velocity(SN['spectra'], lines_dict)
    SN['expansion_velocities'] = spectra_velocity.make_velocity_df(SN, lines_dict)

spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SNiPTF14hls['expansion_velocities']], 'absorption')
spectra_velocity.plot_expansion_velocities([SN2018hmx['expansion_velocities'], SN2018aad['expansion_velocities']], 'absorption')

SN2018hmx['expansion_velocities'].to_csv(r'results\sN2018hmx_expansion_velocities.csv')
SNiPTF14hls['expansion_velocities'].to_csv(r'results\SNiPTF14hls_expansion_velocities.csv')
SN2018aad['expansion_velocities'].to_csv(r'results\SN2018aad_expansion_velocities.csv')

plt.show()

