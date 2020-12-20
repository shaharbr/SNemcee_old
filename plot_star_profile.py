import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt


'''
mass: gram
radius: cm
temp: K
rho: gr / cm^3
vel: cm / s

'''

s17 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's17.0.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

s17_k40_r3000 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's17.0_K40_R3000.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

s9 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's9.0.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

m_sol = 1.989 * (10 ** 33)  # gram
r_sol = 6.9634 * (10 ** 10)  # cm


# print
K = 40
r_vec = s17_k40_r3000['radius']

csm_rho = np.ones(len(r_vec)) * K *(10**17) / (r_vec**2)
r_start = 900 * r_sol
r_end = 3000 * r_sol
csm_rho_integral = 4 * np.pi * K*(10**17) * (r_end - r_start) # volumetric spheric integral
print(csm_rho_integral / m_sol)
csm_rho_integral = 4 * np.pi * K*(10**17) * (r_end - 0) # volumetric spheric integral
print(csm_rho_integral / m_sol)
exit()


# plt.figure()
# plt.plot(s17['radius'], s17['mass'])

# plt.figure()
# plt.plot(s17['radius'], s17['rho'])

# plt.figure()
# plt.plot(np.log10(s17['radius']), s17['rho'])

plt.figure()
plt.plot(s17['mass'] / m_sol, np.log10(s17['rho']))
# plt.plot(s9['mass'] / m_sol, np.log10(s9['rho']))
# plt.plot(s17_k40_r3000['mass'] / m_sol, np.log10(s17_k40_r3000['rho']))
plt.plot(s17_k40_r3000['mass'] / m_sol, np.log10(csm_rho))

plt.figure()
plt.plot(s17['radius'] / r_sol, np.log10(s17['rho']))
# plt.plot(s9['radius'] / r_sol, np.log10(s9['rho']))
# plt.plot(s17_k40_r3000['radius'] / r_sol, np.log10(s17_k40_r3000['rho']))
plt.plot(s17_k40_r3000['radius'] / r_sol, np.log10(csm_rho))


plt.figure()
plt.plot(s17['radius'] / r_sol, s17['mass'] / m_sol)


# plt.figure()
# plt.plot(s17['radius'], s17['vel'])




plt.show()





