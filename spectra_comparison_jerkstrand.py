import pandas as pd
import os
from matplotlib import pyplot as plt

smooth_d357 = pd.read_csv(os.path.join('results','calibrated_d357_spectrum.csv'))
# smooth_d357['index'].rename('plpl')
models_dir = os.path.join('data', 'jerkstrand_models')
# filelist = os.listdir(models_dir)
# filelist.remove('INFO.txt')
filelist = ['mzams12_400d.dat', 'mzams15_369d.dat', 'mzams19_369d.dat']
filenames = ['mzams12 400d', 'mzams15 369d', 'mzams19 369d']


# filelist = ['mzams12_306d.dat']
plt.figure(figsize=(16, 12))

for i in range(len(filelist)):
    model = pd.read_csv(os.path.join(models_dir, filelist[i]), sep=r'\s+', names=['wave', 'flux'])
    plt.plot(model['wave'], model['flux'], label='Jerkstrand+2014 model '+filenames[i])
plt.title('SN 2018hmx 357d vs. Jerkstrand+2014 models', fontsize=16)
plt.plot(smooth_d357['wave'], smooth_d357['flux'], color='k', label='SN2018hmx 357d')
plt.ylim(0, 10e-15)
plt.xlim(4300, 9400)
plt.legend(fontsize=12)
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux (erg s^-1 cm^-2 A^-1)')
plt.tight_layout()
plt.savefig(os.path.join('figures', 'SN2018hmx_Jerkstrand14_models_comparison.png'))

plt.show()