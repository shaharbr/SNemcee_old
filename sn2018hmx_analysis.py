from matplotlib import pyplot as plt
import data_import
import lightcurve_param_fit
import lightcurve
import pandas as pd
import numpy as np
'''
# TODO - ask Iair:
- what should I do regarding the G and o filters which dont have galactic extinction values in the NES?
- removing R filter from bolometric - griffin's code limits to epoches with at least 3 filters by deafult.
If removing R, much of the radioactive tail is excluded, or changing to to at least 2 filters (then high error) 

# TODO additional analysis after GSP: 
- add comparison of lightcurve params to valenti 2D plots
- compare velocities to Gutierrez sample. (and light curves?)
  (https://iopscience.iop.org/article/10.3847/1538-4357/aa8f52,
  in particular Figures 11, and 18, and Table 8
  values from the absorption method not the FWHM one).
- Derive explosion time by comparing temperature and/or B-V colors to other SNe


- check ATLAS and ZTF might have more photometry points on their webpage.

'''


# TODO add lightcurve compared to 18ad?
# TODO add atlas photometry data

plt.rcParams['font.sans-serif'] = 'Arial'

# input correction parameters for all investigated SNe
distance_modulus = {'SN2018hmx': 36.06,
                    'SN1999em': 30.03,
                    'SN2004ej': 32.3,
                    'SN2012A': 29.73 ,
                    'ASASSN14kg': 34.04}

z = {'SN2018hmx': 0.037,
     'SN1999em': 0.0024,
     'SN2004ej': 0.0091,
     'SN2012A': 0.002291429,  # note: SN2012A z calculated as average of multiple z values onon the sn2 catalog
     'ASASSN14kg': 0.014477}

discovery_date = {'SN2018hmx': '2018-10-17 15:30:14',
                  'SN1999em': '1999-10-29 10:33:00',
                  'SN2004ej': '2004-09-10',
                  'SN2012A': '2012-01-07',
                  'ASASSN14kg': '2014-11-17'}

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


# import photometry data files
lco_phot = data_import.lco_phot('sn2018hmx_20181101-20190505_lcophot.txt')
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

def make_alllightcurve_df(SN_dict):
    sources = SN_dict['lightcurve'].keys()
    all_df = []
    for source in sources:
        df = SN_dict['lightcurve'][source]['df'].copy(deep=True)
        df = df.loc[df['mag'] > 0]
        df['source'] = source
        all_df.append(df)
    all_df = pd.concat(all_df, sort=False)
    all_df = all_df[['mjd', 'mag', 'dmag', 'filter', 'source']]
    all_df.reset_index(inplace=True, drop=True)
    return all_df


def save_ascii(dataframe, filename):
    df = dataframe.copy(deep=True)
    df.rename(columns={'mjd': 'MJD', 'filter': 'filt'}, inplace=True)
    df['nondet |'] = 'FALSE |'
    header_row = pd.DataFrame(np.array(df.keys()).reshape(1, 6), columns=list(df.keys()))
    df = pd.concat([header_row, df], ignore_index=True)
    def equalspace(x):
        x = str(x)
        len_x = len(x)
        x = '|' + ' '*(20-len_x) + x
        return x
    for column in df.keys():
        df[column] = df[column].apply(equalspace)
    df.to_csv(filename, sep=str(' '), index=False, header=False)
    return df



SN2018hmx_lightcurves = make_alllightcurve_df(SN2018hmx)
ascii18hmx = save_ascii(SN2018hmx_lightcurves, 'ascii18hmx.ascii')


lightcurve.lightcurve_plot([SN2018hmx, SN1999em], main_SN='SN2018hmx')
lightcurve.lightcurve_plot([SN2018hmx, SN2004ej], main_SN='SN2018hmx')
lightcurve.lightcurve_plot([SN2018hmx, SN2012A], main_SN='SN2018hmx')
lightcurve.lightcurve_plot([SN2018hmx, ASASSN14kg], main_SN='SN2018hmx')

# light curve parameters:
# TODO also for bolometric, not only for the V

# calculate v magnitude at day 50 as the closest observation to day 50
LCO_mag_df = SN2018hmx['lightcurve']['LCO']['df']
vmag_df = LCO_mag_df.loc[LCO_mag_df['filter'] == 'V']
vmag_50 = vmag_df.loc[vmag_df['t_from_discovery']<50, 'abs_mag'].iloc[-1]
vmag_50_err = vmag_df.loc[vmag_df['t_from_discovery']<50, 'dmag'].iloc[-1]


s50V, s50V_sigma = lightcurve_param_fit.calc_s50V(SN2018hmx)
p0, p0_sigma = lightcurve_param_fit.calc_p0(SN2018hmx, time_range=[140, 190])
Vmax, Vmax_err = lightcurve_param_fit.calc_Vmax(SN2018hmx)

lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 'p0')
lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 's50V')

sampler = lightcurve_param_fit.SN_lightcurve_params(SN2018hmx)
lightcurve_param_fit.chain_plots(sampler.chain)
lightcurve_param_fit.plot_v_lightcurve_with_fit(SN2018hmx, sampler)


# TODO what is the correct method for calculating parameter uncertainties? in particular p0 and s50v, but also MCMC.
param_results = lightcurve_param_fit.get_param_results_dict(sampler)
param_results['vmag_50'] = vmag_50
param_results['vmag_50_err'] = vmag_50_err
param_results['s50V'] = s50V
param_results['s50V_err'] = s50V_sigma
param_results['p0'] = p0
param_results['p0_err'] = p0_sigma
param_results['Vmax'] = Vmax
param_results['e_Vmax'] = Vmax_err
param_results['Name'] = 'SN2018hmx'
print(param_results)
param_results = pd.DataFrame(param_results, index=[0])
param_results.to_csv('SN2018hmx_param_results.csv')


lightcurve.BV_plot(SN2018hmx)




plt.show()

