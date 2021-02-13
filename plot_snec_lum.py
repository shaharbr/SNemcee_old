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

plt.figure()
# for source in ['mesa', 'kepler']:
for M in ['19.0']:
    for E in ['1.25']:
        for Mix in ['15.0']:
            for Ni in [0.056]:
                lum = pd.read_csv(os.path.join('data', 'example_lum_snec',
                                             'M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix)+'_R0_K0','lum_observed.dat'),
                                names=['time', 'Lum'], sep=r'\s+')
                lum['time'] = lum['time'] / 86400
                plt.plot(lum['time'], np.log10(lum['Lum']), label='snec model M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix))
                plt.xlim(-2, 180)
                plt.yticks(np.arange(41, 43.5, 0.5))
                plt.xticks(np.arange(0, 180, 50))
                plt.ylim(41, 43)
                plt.xlabel('time (days)')
                plt.ylabel('lum (erg/s, log)')
                # plt.yscale('log')
                # plt.xscale('log')
                # plt.ylim(float(5*10**40), float(3.5*10**42))




sn1999em = pd.read_csv(os.path.join('results', 'bersten_1999em_bol.csv'), names=['time', 'lum'])

plt.errorbar(sn1999em['time'], sn1999em['lum'], marker='o', color='k', linestyle='None', label='SN1999em')
plt.legend()
plt.title('SN1999em vs. SNEC M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix))
plt.tight_layout()

plt.savefig(os.path.join('figures','lum_logY_M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix)+'_R0_K0.png'))

plt.show()
