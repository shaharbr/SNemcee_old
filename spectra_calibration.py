from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
os.environ['PYSYN_CDBS'] = os.path.join('..', 'cdbs')
import pysynphot as S

phot_data = pd.read_csv(os.path.join('results', 'SN2018hmx_Arizona_LCO_mag_df'))

def get_phot_on_day(filter, spect_day):
    filter_phot = phot_data.loc[phot_data['filter'] == filter]
    discovery_day = 58408.646
    z = 0.038
    day_notrest = spect_day - discovery_day
    rest_day = day_notrest / (1 + z)
    print('restday', str(rest_day))
    phot_on_day = np.interp(rest_day, filter_phot['t_from_discovery'], filter_phot['abs_mag'])
    return phot_on_day


val = get_phot_on_day('g', 58494.40432000)
print('phot_on_day', val)

# day
spectra_path = os.path.join('data', 'uncalibrated_18hmx_spectra', '2018hmx_20190111_2458494.90432_1.ascii')



spectra = pd.read_csv(spectra_path, sep=r'\s+', names=["x", "y"])
spectra = spectra.loc[spectra['x'] > 3800]


print(spectra)
plt.figure()
plt.plot(spectra['x'], spectra['y'])

gp_path = os.path.join('data', 'LasCumbres_LasCumbres.SDSS_gp.dat')
gp = pd.read_csv(gp_path, sep=r'\s+', names=["x", "y"])
plt.figure()
plt.plot(gp['x'], gp['y'])
bp = S.FileBandpass(gp_path)
print(bp.equivwidth())
plt.show()

sp = S.FileSpectrum(spectra_path)
# TODO understand which one is the right mag, if it needs corrections
obs = S.Observation(sp, bp)
print(obs.countrate())
print('stmag', obs.effstim('stmag'))
print('abmag', obs.effstim('abmag'))
print('obmag', obs.effstim('obmag'))
print('vegamag', obs.effstim('vegamag'))

print(obs.efflam())
print(obs.efflam(binned=False))

plt.show()
exit()
