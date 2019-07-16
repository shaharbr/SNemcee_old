from matplotlib import pyplot as plt
import data_import
import lightcurve_param_fit
import lightcurve

'''
# TODO:
- standardize the SN_dict object so it is the same for both the photometry and spectroscopy analyses
 TODO turn the data input dicts to stardardized files that can be imported with the data (csv?)


# TODO - ask Iair:
- what should I do regarding the G and o filters which dont have galactic extinction values in the NES?
- Derive explosion time by comparing temperature and/or B-V colors to other SNe


# TODO additional analysis after GSP: 
- add comparison of lightcurve params to valenti 2D plots
- analyse and compare velocities to 18aad (data on SNEx
- compare velocities to Gutierrez sample. (and light curves?)
  (https://iopscience.iop.org/article/10.3847/1538-4357/aa8f52,
  in particular Figures 11, and 18, and Table 7 (and table 8?)
  values from the absorption method not the FWHM one).

- exclude  the r-band data in the blackbody fits,
- Derive explosion time by extrapolating blackbody radius back in time to 0 radius.
- check ATLAS might have more photometry points on their webpage.

'''




plt.rcParams['font.sans-serif'] = 'Arial'

# input correction parameters for all investigated SNe
distance_modulus = {'SN2018hmx': 36.06, 'SN1999em': 30.03, 'SN2004ej': 32.3, 'SN2012A': 29.73 , 'ASASSN14kg': 34.04}
z = {'SN2018hmx': 0.037, 'SN1999em': 0.0024, 'SN2004ej': 0.0091, 'SN2012A': 0.002291429, 'ASASSN14kg': 0.014477}
# note: SN2012A z calculated as average of multiple z values onon the sn2 catalog
discovery_date = {'SN2018hmx': '2018-10-17 15:30:14', 'SN1999em': '1999-10-29 10:33:00',
                  'SN2004ej': '2004-09-10', 'SN2012A': '2012-01-07', 'ASASSN14kg': '2014-11-17'}
galactic_extinction = {'SN2018hmx': {'U': 0.206, 'B': 0.173, 'V': 0.131, 'R': 0.103, 'I': 0.072,
                       'u': 0.202, 'g': 0.157, 'r': 0.109, 'i': 0.081, 'z': 0.060,
                       'J': 0.034, 'H': 0.021, 'K': 0.014, 'L': 0.007,
                       'G': 0, 'o': 0},
                       'SN1999em': {'U': 0.176, 'B': 0.147, 'V': 0.111, 'R': 0.088, 'I': 0.061,
                       'u': 0.172, 'g': 0.134, 'r': 0.093, 'i': 0.069, 'z': 0.051,
                       'J': 0.029, 'H': 0.018, 'K': 0.012, 'L': 0.006,
                       'G': 0, 'o': 0},
                       'SN2004ej': {'V': 0},
                       'SN2012A': {'B': 0, 'V': 0, 'U': 0, 'UVW1': 0, 'UVW2': 0, 'UVM2': 0},
                       'ASASSN14kg': {'B': 0, 'V': 0, 'U': 0, 'i': 0, 'r': 0, 'g': 0}}
# define colormap for plotting, the colors each filter will be presented in
colormap = {'i': 'firebrick', 'r': 'tomato', 'g': 'turquoise',
            'V': 'limegreen', 'B': 'blue', 'U': 'darkorchid', 'G': 'teal', 'R': 'tomato', 'I': 'firebrick',
            'o': 'orange',
            'UVW1': 'darkorchid', 'UVW2': 'darkorchid', 'UVM2': 'darkorchid'}

# convert discovery dates to MJD
for SN in discovery_date.keys():
    discovery_date[SN] = data_import.convert_to_mjd(discovery_date[SN], from_datetime=True)


# import photometry data files
lco_phot = data_import.lco_phot('lco_photometry.txt')
gaia_phot = data_import.gaia_phot('gaia18doi.txt')
atlas_phot = data_import.atlas_phot('ATLAS1_ACAM1.txt')
ztf_phot_new = data_import.ztf_phot_new('ztf_dr1_lightcurve.txt')
sn1999em_leonard_phot = data_import.leonard_phot('sn1999em_UBVRI_leonard02.txt')
sn2004ej_phot = data_import.sne_catalog_phot('SN2004ej.csv')
sn2012A_phot = data_import.sne_catalog_phot('SN2012A.csv')
ASASSN14kg_phot = data_import.sne_catalog_phot('ASASSN-14kg.csv')




# sort the light curve data of each SN in a dict, which also has the plotting guides for each light curve
lightcurves = {'SN2018hmx': {
                    'Gaia': {'df': gaia_phot, 'marker': 'D', 'Linestyle': 'None'},
                    'ATLAS': {'df': atlas_phot, 'marker': 's', 'Linestyle': 'None'},
                    'ZTF': {'df': ztf_phot_new, 'marker': '^', 'Linestyle': 'None'},
                    'LCO': {'df': lco_phot, 'marker': 'o', 'Linestyle': 'None'}},
                'SN1999em': {
                    'Leonard': {'df': sn1999em_leonard_phot, 'marker': 'None', 'Linestyle': '--'}},
                'SN2004ej': {
                    'Planck': {'df': sn2004ej_phot, 'marker': 'None', 'Linestyle': '--'}},
                'SN2012A': {
                    'Swift': {'df': sn2012A_phot, 'marker': 'None', 'Linestyle': '--'}},
                'ASASSN14kg': {
                    'ASASSN': {'df': ASASSN14kg_phot, 'marker': 'None', 'Linestyle': '--'}},
                }



SN2018hmx = data_import.make_SN_dict('SN2018hmx', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)
SN1999em = data_import.make_SN_dict('SN1999em', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)
SN2004ej = data_import.make_SN_dict('SN2004ej', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)
SN2012A = data_import.make_SN_dict('SN2012A', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)
ASASSN14kg = data_import.make_SN_dict('ASASSN14kg', lightcurves, z, discovery_date, distance_modulus, galactic_extinction)


SN2018hmx = lightcurve.add_rest_frame_days_from_discovery(SN2018hmx)
SN2004ej = lightcurve.add_rest_frame_days_from_discovery(SN2004ej)
SN2012A = lightcurve.add_rest_frame_days_from_discovery(SN2012A)
ASASSN14kg = lightcurve.add_rest_frame_days_from_discovery(ASASSN14kg)


for SN in [SN2018hmx, SN1999em, SN2012A, SN2004ej, ASASSN14kg]:
    SN = lightcurve.remove_glactic_extinction(SN)
    SN = lightcurve.add_absolute_magnitude(SN)

SN2018hmx = lightcurve.remove_data_before_discovery(SN2018hmx)
SN1999em = lightcurve.remove_data_before_discovery(SN1999em)


def remove_LCO_outlier(SN_dict):
    df = SN_dict['lightcurve']['LCO']['df']
    df = df.loc[df['t_from_discovery'] < 190]
    return df

SN2018hmx['lightcurve']['LCO']['df'] = remove_LCO_outlier(SN2018hmx)


lightcurve.lightcurve_plot([SN2018hmx, SN1999em], main_SN='SN2018hmx')
lightcurve.lightcurve_plot([SN2018hmx, SN2004ej], main_SN='SN2018hmx')
lightcurve.lightcurve_plot([SN2018hmx, SN2012A], main_SN='SN2018hmx')
lightcurve.lightcurve_plot([SN2018hmx, ASASSN14kg], main_SN='SN2018hmx')

# light curve parameters:

s50V_2018hmx = lightcurve_param_fit.calc_s50V(SN2018hmx)
p0_2018hmx = lightcurve_param_fit.calc_p0(SN2018hmx, time_range=[140, 190])

# lightcurve_plot([SN2018hmx], main_SN='SN2018hmx', V50_line=v50_regression_params)

# lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 'p0')
# lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 's50V')

# sampler = lightcurve_param_fit.SN_lightcurve_params(SN2018hmx)
# lightcurve_param_fit.chain_plots(sampler.chain)
# lightcurve_param_fit.plot_v_lightcurve_with_fit(SN2018hmx, sampler)


# print(sampler.chain.shape)




plt.show()

