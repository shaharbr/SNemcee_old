from matplotlib import pyplot as plt
import pandas as pd
import os

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

plt.figure()
for M in ['13.0']:
    for E in ['0.9']:
        for Mix in ['2.0', '5.0', '8.0']:

            lum = pd.read_csv(os.path.join('data', 'example_lum_snec',
                                         'M'+str(M)+'_Ni0.07_E'+str(E)+'_Mix'+Mix+'_R0.1_K0.1','lum_observed.dat'),
                            names=['time', 'Lum'], sep=r'\s+')
            # energies = energies.iloc[1:]
            # convert sec to days
            lum['time'] = lum['time'] / 86400
            plt.plot(lum['time'], lum['Lum'], label=Mix)
            plt.legend()
            plt.xlim(-2, 200)
            plt.xlabel('time (days, log)')
            plt.ylabel('lum (erg/s, log)')
            plt.yscale('log')
            # plt.xscale('log')
            # plt.ylim(float(5*10**40), float(3.5*10**42))
            plt.title('M'+str(M)+'_Ni0.07_E'+str(E)+'_Mix'+Mix+'_R0.1_K0.1')
            plt.tight_layout()
            plt.savefig(os.path.join('figures',
                                         'lum_logY_M'+str(M)+'_Ni0.07_E'+str(E)+'_Mix'+Mix+'_R0.1_K0.1.png'))

plt.show()
