


exit()

import pandas as pd
import os
from matplotlib import pyplot as plt
import numpy as np

# 2018hmx
smooth_d357 = pd.read_csv(os.path.join('results','calibrated_d357_spectrum.csv'))
SN18hmx_avg = smooth_d357['flux'].mean()
print(SN18hmx_avg)



# maguire 2014et, normalized to 18hmx
SN2014et = pd.read_csv(os.path.join('data','SN2004et_spect_369d.csv'))
SN2014et.rename({'wavelength': 'wave'}, axis='columns', inplace=True)
z = 0.0011025
SN2014et['wave'] = SN2014et['wave'] / (1+z)
SN2014et['flux'] = pd.Series(SN2014et['flux']).rolling(7, center=True).mean()
SN14et_avg = SN2014et['flux'].mean()
SN2014et['flux'] = SN2014et['flux'] * SN18hmx_avg / SN14et_avg
print(SN14et_avg)



# 1999em, normalized to 18hmx
# discovery mjd 51480
# mjd 51899, day 419, rest day 418.000144
SN1999em_d418 = pd.read_csv(os.path.join('data','1999em_spectra','SN1999em_2000-12-21_00-00-00_Lick-3m_KAST_None.ascii')
                       , sep=r'\s+', names=['wave', 'flux'])
# mjd 51813, day 333, rest day 332.205365
SN1999em_d332 = pd.read_csv(os.path.join('data','1999em_spectra','1999em_2000-09-26_00-00-00_Lick-3m_KAST_SUSPECT.dat')
                       , sep=r'\s+', names=['wave', 'flux'])
z = 0.002392
for SN in [SN1999em_d418, SN1999em_d332]:
    SN['wave'] = SN['wave'] / (1+z)
    # SN['flux'] = pd.Series(SN['flux']).rolling(7, center=True).mean()
    SN_avg = SN['flux'].mean()
    SN['flux'] = SN['flux'] * SN18hmx_avg / SN_avg
    print(SN_avg)



# jerkstrand models
models_dir = os.path.join('data', 'jerkstrand_models')
filelist = ['mzams12_400d.dat', 'mzams15_369d.dat']
filenames = ['mzams12 400d', 'mzams15 369d']
times = [400, 369]
colors = ['tab:blue', 'tab:green']


plt.figure(figsize=(16, 12))


for i in range(len(filelist)):
    model = pd.read_csv(os.path.join(models_dir, filelist[i]), sep=r'\s+', names=['wave', 'flux'])

    # multiply flux by 2.5x to account for higher Ni (0.15 instead of 0.06 assumed in the original model)
    # model['flux'] = model['flux'] * 2.5

    # correction for different time in models
    model['flux'] = model['flux'] * np.exp(-(357 - times[i]) / 111.3)

    # plt.plot(model['wave'], model['flux'], label='Jerkstrand+2014 model '+filenames[i], color=colors[i])
plt.title('SN 2018hmx 357d vs. SN 1999em 332d and 418d', fontsize=16)
plt.plot(SN1999em_d332['wave'], SN1999em_d332['flux'], color='tab:blue', label='SN1999em 332d')
plt.plot(SN1999em_d418['wave'], SN1999em_d418['flux'], color='tab:red', label='SN1999em 418d')
plt.plot(smooth_d357['wave'], smooth_d357['flux'], color='k', label='SN2018hmx 357d')


plt.ylim(-0.2e-15, 8.3e-15)
plt.xlim(4300, 9400)
plt.legend(fontsize=12)
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux (erg s^-1 cm^-2 A^-1)')
plt.tight_layout()
plt.savefig(os.path.join('figures', 'SN2018hmx_SN1999em_comparison.png'))


plt.show()