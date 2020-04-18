from matplotlib import pyplot as plt
import data_import
import lightcurve_param_fit
import lightcurve
import pandas as pd
import numpy as np




# import photometry data files for 2018hmx
atlas_phot_o = data_import.atlas_flux('data\SN2018hmx_flux_raw_wightedavg_orange.txt', 'o')
atlas_phot_c = data_import.atlas_flux('data\SN2018hmx_flux_raw_wightedavg_cyan.txt', 'c')
atlas_phot = pd.concat([atlas_phot_o, atlas_phot_c], ignore_index=True)
lco_phot = data_import.lco_phot(r'data\sn2018hmx_20181101-20190505_lcophot.txt')
gaia_phot = data_import.gaia_phot(r'data\gaia18doi.txt')
ztf_phot_new = data_import.ztf_phot_new(r'data\ztf_dr1_lightcurve.txt')
arizona_phot = pd.read_csv(r'data\SN2018hmx_phot_arizona.csv')


# import photometry data files for other comparison SNe
sn1999em_leonard_phot = data_import.leonard_phot(r'data\sn1999em_UBVRI_leonard02.txt')
ASASSN14kg_phot = data_import.sne_catalog_phot(r'data\ASASSN-14kg.csv')
ksp_phot = data_import.ksp_phot('data\KSP-SN-2016kf.txt')
sn2018aad_phot = data_import.sn18aad_phot(r'data\sn2018aad_20180305-20190325_pyzogy_lcophot.txt')

# import correction params
correction_params = pd.concat([
    pd.read_csv(r'data\SN2018hmx_correction.csv'),
    pd.read_csv(r'data\ASASSN14kg_correction.csv'),
    pd.read_csv(r'data\SN1999em_correction.csv'),
    pd.read_csv(r'data\KSP_correction.csv'),
    pd.read_csv(r'data\SN2018aad_correction.csv')
    ], sort=False, ignore_index=True)
correction_params.set_index('Name', inplace=True)
colormap = pd.read_csv(r'data\colormap.csv')

# print(correction_params)

# sort the light curve data of each SN in a dict, which also has the plotting guides for each light curve
SN2018hmx = {'Name': 'SN2018hmx'}
SN2018hmx['lightcurve'] = {'Gaia': {'df': gaia_phot, 'marker': 'o', 'Linestyle': 'None'},
                          'ATLAS': {'df': atlas_phot, 'marker': 'o', 'Linestyle': 'None'},
                          'P60': {'df': ztf_phot_new, 'marker': 'o', 'Linestyle': 'None'},
                          'Las Cumbres': {'df': lco_phot, 'marker': 'o', 'Linestyle': 'None'},
                          'Arizona': {'df': arizona_phot, 'marker': 'o', 'Linestyle': 'None'}}

SN1999em = {'Name': 'SN1999em'}
SN1999em['lightcurve'] = {'Leonard': {'df': sn1999em_leonard_phot, 'marker': 'None', 'Linestyle': '--'}}

ASASSN14kg = {'Name': 'ASASSN14kg'}
ASASSN14kg['lightcurve'] = {'ASASSN': {'df': ASASSN14kg_phot, 'marker': 'None', 'Linestyle': '--'}}

ksp = {'Name': 'KSP-SN-2016kf'}
ksp['lightcurve'] = {'Las Cumbres': {'df': ksp_phot, 'marker': 'None', 'Linestyle': '--'}}

SN2018aad = {'Name': 'SN2018aad'}
SN2018aad['lightcurve'] = {'Las Cumbres': {'df': sn2018aad_phot, 'marker': '^', 'Linestyle': 'None'}}


# processing
SN2018hmx = lightcurve.add_rest_frame_days_from_discovery(SN2018hmx, correction_params.loc['SN2018hmx'])
SN2018hmx = lightcurve.remove_glactic_extinction(SN2018hmx, correction_params.loc['SN2018hmx'])
SN2018hmx = lightcurve.add_absolute_magnitude(SN2018hmx, correction_params.loc['SN2018hmx'])

ASASSN14kg = lightcurve.add_rest_frame_days_from_discovery(ASASSN14kg, correction_params.loc['ASASSN14kg'])
ASASSN14kg = lightcurve.remove_glactic_extinction(ASASSN14kg, correction_params.loc['ASASSN14kg'])
ASASSN14kg = lightcurve.add_absolute_magnitude(ASASSN14kg, correction_params.loc['ASASSN14kg'])

SN1999em = lightcurve.remove_glactic_extinction(SN1999em, correction_params.loc['SN1999em'])
SN1999em = lightcurve.add_absolute_magnitude(SN1999em, correction_params.loc['SN1999em'])

ksp = lightcurve.add_rest_frame_days_from_discovery(ksp, correction_params.loc['KSP-SN-2016kf'])
ksp = lightcurve.remove_glactic_extinction(ksp, correction_params.loc['KSP-SN-2016kf'])
ksp = lightcurve.add_absolute_magnitude(ksp, correction_params.loc['KSP-SN-2016kf'])

SN2018aad = lightcurve.add_rest_frame_days_from_discovery(SN2018aad, correction_params.loc['SN2018aad'])
SN2018aad = lightcurve.remove_glactic_extinction(SN2018aad, correction_params.loc['SN2018aad'])
SN2018aad = lightcurve.add_absolute_magnitude(SN2018aad, correction_params.loc['SN2018aad'])


SN2018hmx_lightcurves = data_import.make_alllightcurve_df(SN2018hmx)
# ascii18hmx = data_import.save_ascii(SN2018hmx_lightcurves, r'data\ascii18hmx.ascii')

SN2018aad_lightcurves = data_import.make_alllightcurve_df(SN2018aad)
ascii18aad = data_import.save_ascii(SN2018aad_lightcurves, r'data\ascii18aad.ascii')

lightcurve.lightcurve_plot_shift(ksp, correction_params)

# lightcurve.lightcurve_plot([SN2018hmx], 'SN2018hmx', correction_params)
# lightcurve.lightcurve_plot([SN2018aad], 'SN2018aad', correction_params)
# lightcurve.lightcurve_plot_shift(SN2018hmx, correction_params)



lightcurve.lightcurve_plot([SN2018hmx], 'SN2018hmx', correction_params)
lightcurve.lightcurve_plot([SN2018aad], 'SN2018aad', correction_params)
# lightcurve.lightcurve_plot([SN2018hmx, ASASSN14kg], 'SN2018hmx', correction_params)


# light curve parameters:
# TODO also for bolometric, not only for the V

# calculate v magnitude at day 50 as the closest observation to day 50

print('SN2018hmx')
LCO_mag_df = SN2018hmx['lightcurve']['Las Cumbres']['df']
vmag_df = LCO_mag_df.loc[LCO_mag_df['filter'] == 'V']
vmag_50 = vmag_df.loc[vmag_df['t_from_discovery']<50, 'abs_mag'].iloc[-1]
vmag_50_err = vmag_df.loc[vmag_df['t_from_discovery']<50, 'dmag'].iloc[-1]
s50V, s50V_sigma = lightcurve_param_fit.calc_s50V(SN2018hmx)
p0, p0_sigma = lightcurve_param_fit.calc_p0(SN2018hmx, time_range=[140, 190])
Vmax, Vmax_err = lightcurve_param_fit.calc_Vmax(SN2018hmx)
# lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 'p0')
lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018hmx, 's50V')

print(s50V, s50V_sigma)
print(p0, p0_sigma)
print(Vmax, Vmax_err)
print(vmag_50)
print(vmag_50_err)

# TODO the p0 is manually inserted in the lightcurve param fit code, so it can affect fit values
sampler = lightcurve_param_fit.SN_lightcurve_params(SN2018hmx)
lightcurve_param_fit.chain_plots(sampler.chain)
lightcurve_param_fit.plot_v_lightcurve_with_fit(SN2018hmx, sampler)

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
param_results.to_csv('results\SN2018hmx_param_results.csv')

lightcurve.BV_plot(SN2018hmx)



print('SN2018aad')
LCO_mag_df = SN2018aad['lightcurve']['Las Cumbres']['df']
vmag_df = LCO_mag_df.loc[LCO_mag_df['filter'] == 'V']
vmag_50 = vmag_df.loc[vmag_df['t_from_discovery']<50, 'abs_mag'].iloc[-1]
vmag_50_err = vmag_df.loc[vmag_df['t_from_discovery']<50, 'dmag'].iloc[-1]
s50V, s50V_sigma = lightcurve_param_fit.calc_s50V(SN2018aad)
p0, p0_sigma = lightcurve_param_fit.calc_p0(SN2018aad, time_range=[140, 300])
Vmax, Vmax_err = lightcurve_param_fit.calc_Vmax(SN2018aad)
# lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018aad, 'p0')
lightcurve_param_fit.plot_v_lightcurve_with_slope(SN2018aad, 's50V')

print(s50V, s50V_sigma)
print(p0, p0_sigma)
print(Vmax, Vmax_err)
print(vmag_50)
print(vmag_50_err)

sampler = lightcurve_param_fit.SN_lightcurve_params(SN2018aad)
lightcurve_param_fit.chain_plots(sampler.chain)
lightcurve_param_fit.plot_v_lightcurve_with_fit(SN2018aad, sampler)

param_results = lightcurve_param_fit.get_param_results_dict(sampler)
param_results['vmag_50'] = vmag_50
param_results['vmag_50_err'] = vmag_50_err
param_results['s50V'] = s50V
param_results['s50V_err'] = s50V_sigma
param_results['p0'] = p0
param_results['p0_err'] = p0_sigma
param_results['Vmax'] = Vmax
param_results['e_Vmax'] = Vmax_err
param_results['Name'] = 'SN2018aad'
print(param_results)
param_results = pd.DataFrame(param_results, index=[0])
param_results.to_csv('results\SN2018aad_param_results.csv')

lightcurve.BV_plot(SN2018aad)



# print('KSP')
# ksp_mag_df = ksp['lightcurve']['Las Cumbres']['df']
# ksp_vmag_df = ksp_mag_df.loc[ksp_mag_df['filter'] == 'V']
# ksp_vmag_50 = ksp_vmag_df.loc[ksp_vmag_df['t_from_discovery']<50, 'abs_mag'].iloc[-1]
# ksp_vmag_50_err = ksp_vmag_df.loc[ksp_vmag_df['t_from_discovery']<50, 'dmag'].iloc[-1]
# s50V, s50V_sigma = lightcurve_param_fit.calc_s50V(ksp)
# Vmax, Vmax_err = lightcurve_param_fit.calc_Vmax(ksp)
# lightcurve_param_fit.plot_v_lightcurve_with_slope(ksp, 's50V')
# print(ksp_vmag_50)
# print(ksp_vmag_50_err)
# print(s50V, s50V_sigma)
# print(Vmax, Vmax_err)



plt.show()

