from matplotlib import pyplot as plt
import pandas as pd
import os
import numpy as np


# for M in ['9.0', '13.0', '15.0',  '17.0']:
for M in ['15.0']:
    for E in ['1.0']:
        energies = pd.read_csv(os.path.join('data', 'example_energy_transfer_snec',
                                     'M'+M+'_Ni0.02_E'+E+'_Mix1.0_R600_K0.001','conservation.dat'),
                        names=['time','E_gravitational','E_internal','E_kinetic','E_total','EtotmInt'], sep=r'\s+')
        # energies = energies.iloc[1:]
        # convert sec to days
        energies['time'] = energies['time'] / 86400
        energies['E_gravitational'] = -energies['E_gravitational']
        plt.figure()
        for Etype in ['E_gravitational','E_internal','E_kinetic','E_total']:
            plt.plot(energies['time'], energies[Etype], label=Etype)
        plt.legend()
        plt.xlim(0.01, 200)
        plt.xlabel('time (days, log)')
        plt.ylabel('energy (erg, log)')
        # plt.yscale('log')
        # plt.xscale('log')
        plt.title('M'+M+'_Ni0.02_E'+E+'_Mix1.0_R600_K0.001')
        plt.savefig(os.path.join('figures', 'energy_transfer_evolution_logXlogY_M'+M+'_Ni0.02_E'+E+'_Mix1.0_R600_K0.001.png'))

plt.show()
