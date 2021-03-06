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

s15_kepler = pd.read_csv(os.path.join('data', 'example_star_profiles', 's15.0.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

s15_mesa = pd.read_csv(os.path.join('data', 'example_star_profiles', '15Msol_RSG.short'), sep=r'\s+', skiprows=1,
                  names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

# s9_R0_K0 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's9.0_K1e-10_R0.short'), sep=r'\s+', skiprows=1,
#                   names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])

# s9_K100_R0 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's9.0_K100_R0.short'), sep=r'\s+', skiprows=1,
#                   names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])
#
# s9_K0_R2000 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's9.0_K1e-10_R2000.short'), sep=r'\s+', skiprows=1,
#                   names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])
#
# s9_K100_R2000 = pd.read_csv(os.path.join('data', 'example_star_profiles', 's9.0_K100_R2000.short'), sep=r'\s+', skiprows=1,
#                   names=['i', 'mass', 'radius', 'temp', 'rho', 'vel'])
#
m_sol = 1.989 * (10 ** 33)  # gram
r_sol = 6.9634 * (10 ** 10)  # cm

plt.figure()
plt.plot(s15_kepler['radius'] / r_sol, s15_kepler['rho'], color='k')
plt.ylabel('log rho')
plt.xlabel('radius (in R_sol)')
plt.title('s15_kepler')
# plt.ylim(bottom=10**(-16))
plt.yscale('log')

plt.figure()
plt.plot(s15_kepler['mass'] / m_sol, s15_kepler['rho'], color='k')
plt.ylabel('log rho')
plt.xlabel('mass (in M_sol)')
plt.title('s15_kepler')
# plt.ylim(bottom=10**(-16))
plt.yscale('log')


plt.figure()
plt.plot(s15_mesa['radius'] / r_sol, s15_mesa['rho'], color='k')
plt.ylabel('log rho')
plt.xlabel('radius (in R_sol)')
plt.title('s15_mesa')
# plt.ylim(bottom=10**(-16))
plt.yscale('log')

plt.figure()
plt.plot(s15_mesa['mass'] / m_sol, s15_mesa['rho'], color='k')
plt.ylabel('log rho')
plt.xlabel('mass (in M_sol)')
plt.title('s15_mesa')
# plt.ylim(bottom=10**(-16))
plt.yscale('log')


# print
# K = 40
# r_vec = s9_K100_R2000['radius']

# csm_rho = np.ones(len(r_vec)) * K *(10**17) / (r_vec**2)
# r_start = 900 * r_sol
# r_end = 3000 * r_sol
# csm_rho_integral = 4 * np.pi * K*(10**17) * (r_end - r_start) # volumetric spheric integral
# print(csm_rho_integral / m_sol)
# csm_rho_integral = 4 * np.pi * K*(10**17) * (r_end - 0) # volumetric spheric integral
# print(csm_rho_integral / m_sol)


# plt.figure()
# plt.plot(s9['radius'], s9['mass'])

# plt.figure()
# plt.plot(s9['radius'], s9['rho'])

# plt.figure()
# plt.plot(np.log10(s9['radius']), s9['rho'])

# plt.figure()
# plt.plot(s9['mass'] / m_sol, np.log10(s9['rho']), marker='o')
# plt.plot(s9_R0_K0['mass'] / m_sol, np.log10(s9_R0_K0['rho']), marker='o')
# plt.plot(s9_K0_R2000['mass'] / m_sol, np.log10(s9_K0_R2000['rho']), marker='o')
# plt.plot(s9_K100_R0['mass'] / m_sol, np.log10(s9_K100_R0['rho']), marker='o')
# plt.plot(s9_K100_R2000['mass'] / m_sol, np.log10(s9_K100_R2000['rho']), marker='o')
# plt.ylabel('log rho')
# plt.xlabel('mass')


plt.show()





