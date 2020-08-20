from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
os.environ['PYSYN_CDBS'] = os.path.join('..', 'cdbs')
import pysynphot as S
from astropy.io import ascii

filt_colors = dict(zip(['B', 'g', 'V', 'r', 'i'], ['#1700ff', '#00bbff', '#77ff00', '#ff6b00', '#9b0000']))

# filter bandpass functions
gp_path = os.path.join('data', 'LasCumbres_LasCumbres.SDSS_gp.dat')
rp_path = os.path.join('data', 'LasCumbres_LasCumbres.SDSS_rp.dat')
ip_path = os.path.join('data', 'LasCumbres_LasCumbres.SDSS_ip.dat')
B_path = os.path.join('data', 'LasCumbres_LasCumbres.Bessel_B.dat')
V_path = os.path.join('data', 'LasCumbres_LasCumbres.Bessel_V.dat')

# correction parameters for the SN analyzed
correction_table = pd.read_csv(os.path.join('data', 'SN2018hmx_correction.csv'))
B_extinc = float((correction_table['B']))
V_extinc = float((correction_table['V']))
BV_extinc = B_extinc - V_extinc
discovery_mjd = float(correction_table['discovery_date'])
z = float(correction_table['z'])

# spectra data
spectra_d82_path = os.path.join('data', 'uncalibrated_18hmx_spectra', '2018hmx_20190111_2458494.90432_1.ascii')
spectra_d82_mjd = 58494.90432
spectra_d357_path = os.path.join('data', 'uncalibrated_18hmx_spectra', '2018hmx_20191023_2458780.09736_1.ascii')
spectra_d357_mjd = 58780.09736

# photometry data
phot_data = pd.read_csv(os.path.join('results', 'SN2018hmx_Arizona_LCO_mag_df'))

# plt.figure()
# filterlist = phot_data['filter'].unique()
# for filter in filterlist:
#     filt_df = phot_data.loc[phot_data['filter'] == filter]
#     plt.scatter(filt_df['t_from_discovery'], filt_df['abs_mag'], label=filter)
#     plt.legend()
#     plt.gca().invert_yaxis()



# load uncalibrated spectrum from ascii file, and turn into standard spectrum object
def load_spectrum(path, BV_extinction, wavelength_cutoff=False, day=False):
    tab = ascii.read(path, names=['wave', 'flux'])
    if wavelength_cutoff:
        wave = tab['wave'][tab['wave'] > wavelength_cutoff]
        flux = tab['flux'][tab['wave'] > wavelength_cutoff]
    else:
        wave = tab['wave']
        flux = tab['flux']
    sp = S.ArraySpectrum(wave=wave, flux=flux, waveunits='angstrom', fluxunits='flam', keepneg=False)
    sp_ext = sp * S.Extinction(BV_extinction, 'lmcavg')
    if day:
        sp_ext.name = day
    plt.figure()
    plt.plot(sp.wave, sp.flux, 'b', label='E(B-V)=0')
    plt.plot(sp_ext.wave, sp_ext.flux, 'r', label='E(B-V) calc')
    plt.xlabel(sp.waveunits)
    plt.ylabel(sp.fluxunits)
    plt.title(sp.name)
    plt.legend(loc='best')
    return sp_ext

# load bandpass functions from ascii file, and turn into standard bandpass object
def load_bandpass(path, filter):
    tab = ascii.read(path, names=['wave', 'throughput'])
    wave = tab['wave']
    throughput = tab['throughput']
    bandpass = S.ArrayBandpass(wave=wave, throughput=throughput, waveunits='angstrom', name=filter)
    plt.plot(bandpass.wave, bandpass.throughput, label=bandpass.name, color=filt_colors[filter])
    plt.xlabel(bandpass.waveunits)
    plt.ylabel('throughput')
    plt.title(bandpass.name)
    return bandpass

def get_rest_day(spectrum_mjd, discovery_mjd, z):
    day_notrest = spectrum_mjd - discovery_mjd
    rest_day = day_notrest / (1 + z)
    return rest_day

def get_phot_on_day(phot_data, filter, rest_day):
    filter_phot = phot_data.loc[phot_data['filter'] == filter]
    phot_on_day = np.interp(rest_day, filter_phot['t_from_discovery'], filter_phot['mag'])
    return phot_on_day

def get_effstim(spectrum, bandpass):
    obs = S.Observation(spectrum, bandpass, force='taper')
    if bandpass.name in ['g', 'r', 'i']: # small letter filters
        mag_units = 'abmag'
    elif bandpass.name in ['B', 'V']:# capital letter filters
        mag_units = 'vegamag'
    else:
        raise('bandpass can only be one of: g, r, i, B, V')
    effstim = obs.effstim(mag_units)
    return effstim


def get_filter_center(bandpass):
    xm = np.sum(bandpass.wave * bandpass.throughput)
    m = np.sum(bandpass.throughput)
    center = xm / m
    return center


def get_ratio_table(phot_data, spectrum, bandpass_list):
    df = {'ratio': [], 'wave': []}
    for bp in bandpass_list:
        real_phot = get_phot_on_day(phot_data, bp.name, spectrum.name)
        spectra_phot = get_effstim(spectrum, bp)
        ratio = real_phot / spectra_phot
        filt_center = get_filter_center(bp)
        df['ratio'].append(ratio)
        df['wave'].append(filt_center)
    df = pd.DataFrame(df).sort_values('wave').reset_index()
    return df

def fit_libration_line(phot_data, spectrum, bandpass_list):
    df = get_ratio_table(phot_data, spectrum, bandpass_list)
    regresison_params = np.polyfit(df['wave'], df['ratio'], deg=1)
    return regresison_params

d82_day = get_rest_day(spectra_d82_mjd, discovery_mjd, z)
d357_day = get_rest_day(spectra_d357_mjd, discovery_mjd, z)

spectra_d82 = load_spectrum(spectra_d82_path, BV_extinc, 3800, d82_day)
spectra_d357 = load_spectrum(spectra_d357_path, BV_extinc, 4600, d357_day)

plt.figure()
g_bandpass = load_bandpass(gp_path, 'g')
r_bandpass = load_bandpass(rp_path, 'r')
i_bandpass = load_bandpass(ip_path, 'i')
B_bandpass = load_bandpass(B_path, 'B')
V_bandpass = load_bandpass(V_path, 'V')
plt.legend(loc='best')
plt.title('Spectra calibration filters')
plt.savefig(os.path.join('figures', 'spectra_calibration_filters'))


df82 = get_ratio_table(phot_data, spectra_d82, [g_bandpass, r_bandpass, i_bandpass, B_bandpass, V_bandpass])
line82 = fit_libration_line(phot_data, spectra_d82, [g_bandpass, r_bandpass, i_bandpass, B_bandpass, V_bandpass])

df357 = get_ratio_table(phot_data, spectra_d357, [g_bandpass, r_bandpass, i_bandpass, B_bandpass, V_bandpass])
line357 = fit_libration_line(phot_data, spectra_d357, [g_bandpass, r_bandpass, i_bandpass, B_bandpass, V_bandpass])



list_filters = ['B', 'g', 'V', 'r', 'i']
colors = ['#1700ff', '#00bbff', '#77ff00', '#ff6b00', '#9b0000']

plt.figure()
for i in range(5):
    plt.scatter(df82['wave'][i], df82['ratio'][i], label=list_filters[i], color=colors[i])
plt.plot(df82['wave'], line82[1] + line82[0] * df82['wave'], label='linear fit', color='k', alpha=0.5)
plt.title('ratio between the real photometry mag\nand the photometry calculated from the spectrum\n'+
          'on day 82')
plt.xlabel('wavelength (A)')
plt.ylabel('real phot / spectra phot')
plt.legend()
plt.ylim(0.8, 1.4)
plt.tight_layout()
plt.savefig(os.path.join('figures', 'spectra_calibration_ratio_d82'))


plt.figure()
for i in range(5):
    plt.scatter(df357['wave'][i], df357['ratio'][i], label=list_filters[i], color=colors[i])
plt.plot(df357['wave'], line357[1] + line357[0] * df357['wave'], label='linear fit', color='k', alpha=0.5)
plt.title('ratio between the real photometry mag\nand the photometry calculated from the spectrum\n'+
          'on day 357')
plt.xlabel('wavelength (A)')
plt.ylabel('real phot / spectra phot')
plt.legend()
plt.ylim(0.8, 1.4)
plt.tight_layout()
plt.savefig(os.path.join('figures', 'spectra_calibration_ratio_d357'))

calibrated_d82 = spectra_d82.flux * (line82[1] + line82[0] * spectra_d82.wave)
calibrated_d357 = spectra_d357.flux * (line357[1] + line357[0] * spectra_d357.wave)

plt.figure(figsize=(10, 5))
plt.plot(spectra_d82.wave, spectra_d82.flux, label='before_calibration')
plt.plot(spectra_d82.wave, calibrated_d82, label='after_calibration')
plt.xlabel('wavelength (A)')
plt.ylabel('flux (flam)')
plt.title('Day 82 spectrum before and after calibration')
plt.legend()
plt.savefig(os.path.join('figures', 'spectra_calibration_beforeafter_d82'))

plt.figure(figsize=(10, 5))
plt.plot(spectra_d357.wave, spectra_d357.flux, label='before_calibration')
plt.plot(spectra_d357.wave, calibrated_d357, label='after_calibration')
plt.xlabel('wavelength (A)')
plt.ylabel('flux (flam)')
plt.title('Day 357 spectrum before and after calibration')
plt.legend()
plt.savefig(os.path.join('figures', 'spectra_calibration_beforeafter_d357'))


smooth_d82 = pd.Series(spectra_d82.flux).rolling(10, center=True).mean()
smooth_calib_d82 = pd.Series(calibrated_d82).rolling(10, center=True).mean()
plt.figure(figsize=(10, 5))
plt.plot(spectra_d82.wave, smooth_d82, label='before_calibration')
plt.plot(spectra_d82.wave, smooth_calib_d82, label='after_calibration')
plt.xlabel('wavelength (A)')
plt.ylabel('flux (flam)')
plt.title('Day 82 spectrum before and after calibration (smoothed)')
plt.legend()
plt.savefig(os.path.join('figures', 'spectra_calibration_beforeafter_d82_smooth'))

smooth_d357 = pd.Series(spectra_d357.flux).rolling(40, center=True).mean()
smooth_calib_d357 = pd.Series(calibrated_d357).rolling(40, center=True).mean()
plt.figure(figsize=(10, 5))
plt.plot(spectra_d357.wave, smooth_d357, label='before_calibration')
plt.plot(spectra_d357.wave, smooth_calib_d357, label='after_calibration')
plt.xlabel('wavelength (A)')
plt.ylabel('flux (flam)')
plt.title('Day 357 spectrum before and after calibration (smoothed)')
plt.legend()
plt.savefig(os.path.join('figures', 'spectra_calibration_beforeafter_d357_smooth'))



# TODO filters do not have a defined binset in the wavecat table. The waveset of the spectrum will be used instead.

# TODO send iair normalizaiton graphs, and calibrated spectra. ask about units, and about
plt.show()
