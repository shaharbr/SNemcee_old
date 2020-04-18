from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


# Open the given fits file
z = 0.038

filenames = ['2018hmx_20190111_2458494.90432_1', '2018hmx_20191023_2458780.09736_1']
dates = [2458494.90432, 2458780.09736]
dates = [i - 2400000.5 for i in dates]
t_from_discovery = [(i - 58408.646)/(1 + z) for i in dates]

print(t_from_discovery)


def smooth_spectrum(spectrum_y):
    var_list = list(spectrum_y)
    var_list = [x for x in var_list if str(x) != 'nan']
    var_list = var_list[0:100]
    start_variability = np.std(var_list) / 2
    smoothing_window = int(len(spectrum_y)/3000 * start_variability)+1
    smooth_y = spectrum_y.rolling(smoothing_window, center=True).mean()
    return smooth_y

fig, ax = plt.subplots(1, figsize=(9, 6))

for i in range(2):
    hdulist = fits.open('data/'+filenames[i]+'.fits')
    # print(hdulist[0].data*10000000000000000000000000000000)

    scidata = hdulist[0].data
    scidata = np.transpose(scidata)
    scidata = pd.DataFrame(scidata)
    print(scidata)
    # save your new file as a csv
    # np.savetxt('data/'+filenames[i]+'.csv', scidata, delimiter=',')
    avg = np.mean(scidata)
    scidata = scidata / avg - i * 6
    scidata[0] = smooth_spectrum(scidata[0])
    ax.plot(scidata, label='day '+str(int(t_from_discovery[i])))
    ax.set_title('Keck spectra of SN 2018hmx', fontsize=18)
    ax.set_ylabel('Normalized fλ + shift', fontsize=16)


ax.legend(fontsize=14)

plt.show()