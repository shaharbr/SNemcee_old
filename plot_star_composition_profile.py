import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt

first_line = ['ncells','nisotopes']
second_line = ['1.0d0 1.0d0 4.0d0 12.0d0 16.0d0 20.0d0 24.0d0 28.0d0 32.0d0 36.0d0 40.0d0 44.0d0 48.0d0 52.0d0 56.0d0 ']
second_line_elements = ['H_p', 'H', 'He', 'C', 'O', 'Ne', 'Mg', 'Si', 'S', 'Cl', 'Ca', 'Sc', 'Ti', 'Cr', 'Fe']
third_line = ["0.0d0 1.0d0 2.0d0 6.0d0 8.0d0 10.0d0 12.0d0 14.0d0 16.0d0 18.0d0 20.0d0 22.0d0 24.0d0 26.0d0 28.0d0 "]
fourth_line = ["cell mass","cell radius"]+["fraction_"+second_line_elements[i] for i in range(len(second_line_elements))]

s19_path = os.path.join('data', 'example_star_profiles', 's19.0.iso.dat')
s19 = pd.read_csv(s19_path, skiprows=3, sep=r'\s+', names=fourth_line)

R_sol = 6.957* (10**10) # in cm
M_sol = 1.98847 * (10**33)  # in gram
Ni_frac = 0.056 / 15


plt.figure()
plt.xticks([0, 5, 10, 15])
plt.minorticks_on()
plt.xlim(0, 19)
plt.grid(which='both', alpha=0.2)

plt.plot(s19['cell mass'] / M_sol, s19['fraction_H'], label='H', color='red')
plt.plot(s19['cell mass'] / M_sol, s19['fraction_He'], label='He', color='blue')
plt.plot(s19['cell mass'] / M_sol, s19['fraction_C'], label='C', color='purple')
plt.plot(s19['cell mass'] / M_sol, s19['fraction_O'], label='O', color='pink')
plt.plot([0, 15, 15.00001, 16], [20 * Ni_frac, 20 * Ni_frac, 0, 0], label='Ni x 20', color='cyan')
plt.ylabel('mass fraction')
plt.xlabel('Mass (in M_sol)')
plt.legend()
# plt.xscale('log')
# plt.yscale('log')

plt.show()