import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt

s17 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's17.0_K120_R600.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])
print(s17)
csm_rho = np.ones(len(s17['radius'])) * 120 *(10**17) / (s17['radius']**2)

# plt.figure()
# plt.plot(s17['radius'], s17['mass'])

plt.figure()
plt.plot(s17['radius'], s17['rho'])

# plt.figure()
# plt.plot(np.log10(s17['radius']), s17['rho'])

plt.figure()
plt.plot(s17['radius'], np.log10(s17['rho']))
plt.plot(s17['radius'], np.log10(csm_rho))

# plt.figure()
# plt.plot(s17['radius'], s17['vel'])




plt.show()





