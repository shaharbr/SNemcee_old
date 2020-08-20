import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

filepath = os.path.join('..', 'all_snec_data', 'M13.0_Ni0.1_E1.5_Mix3.0_R600_K0.001', 'magnitudes.dat')

mag = pd.read_csv(filepath, names=['time', 'Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'], sep=r'\s+')
# convert seconds to days
mag['time'] = mag['time'] / 86400

fig, ax = plt.subplots()
for filt in ['PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
    ax.plot(mag['time'], mag[filt], label=filt)
    ax.invert_yaxis()
    ax.legend()


plt.show()

print(mag)