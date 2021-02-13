from matplotlib import pyplot as plt
import pandas as pd
import os
import numpy as np

# plt.figure()
# for Mix in ['2.0', '8.0']:
#     lum = pd.read_csv(os.path.join('data', 'example_lum_snec',
#                                  'M9.0_Ni0.17_E0.7_Mix'+Mix+'_R3000_K120','lum_observed.dat'),
#                     names=['time', 'Lum'], sep=r'\s+')
#     energies = energies.iloc[1:]
#     convert sec to days
    # lum['time'] = lum['time'] / 86400
    # plt.plot(lum['time'], lum['Lum'], label='Mix '+str(Mix))
# plt.legend()
# plt.xlabel('time (days)')
# plt.ylabel('lum (erg/s)')
# plt.ylim(float(5*10**40), float(3.5*10**42))
# plt.title('M9.0_Ni0.17_E0.7_R3000_K120')
# plt.savefig(os.path.join('figures',
#                              'lum_M9.0_Ni0.17_E0.7_R3000_K120.png'))

m_sol = 1.989 * (10 ** 33)  # gram

for source in ['mesa', 'kepler']:
    for M in ['15.0']:
        for E in ['1.0']:
            for Mix in ['5.0']:
                for Ni in [0.05]:
                    plt.figure()
                    c_dir = os.path.join('data', 'example_star_composition',
                                                 'M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix)+'_R0_K0_'+source)
                    grid_mass = pd.read_csv(os.path.join(c_dir, 'delta_mass_initial.dat'), names=['gridcell', 'delta_mass'], sep=r'\s+')
                    H_frac = pd.read_csv(os.path.join(c_dir, 'H_init_frac.dat'), names=['gridcell', 'fraction'], sep=r'\s+')
                    He_frac = pd.read_csv(os.path.join(c_dir, 'He_init_frac.dat'), names=['gridcell', 'fraction'], sep=r'\s+')
                    O_frac = pd.read_csv(os.path.join(c_dir, 'O_init_frac.dat'), names=['gridcell', 'fraction'], sep=r'\s+')
                    C_frac = pd.read_csv(os.path.join(c_dir, 'C_init_frac.dat'), names=['gridcell', 'fraction'], sep=r'\s+')
                    Ni_frac = pd.read_csv(os.path.join(c_dir, 'Ni_init_frac.dat'), names=['gridcell', 'fraction'], sep=r'\s+')
                    grid_mass['mass'] = grid_mass['delta_mass'].cumsum() / m_sol
                    plt.plot(grid_mass['mass'], H_frac['fraction'], label='H', color='blue')
                    plt.plot(grid_mass['mass'], He_frac['fraction'], label='He', color='red')
                    plt.plot(grid_mass['mass'], C_frac['fraction'], label='C', color='green')
                    plt.plot(grid_mass['mass'], O_frac['fraction'], label='O', color='black')
                    plt.plot(grid_mass['mass'], Ni_frac['fraction'], label='Ni', color='cyan')
                    plt.legend()
                    plt.yticks(np.arange(0.0, 1.0, 0.2))
                    plt.ylim(0.0, 1.0)
                    plt.xlabel('Mass (in M_sol)')
                    plt.ylabel('Mass fraction')
                    plt.title('M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix)+'_'+source)
                    plt.tight_layout()
                    plt.savefig(os.path.join('figures',
                                             'star_composition_M' + str(M) + '_Ni' + str(Ni) + '_E' + str(E) + '_Mix' + str(
                                                 Mix) + '_R0_K0_'+source+'.png'))

plt.show()
