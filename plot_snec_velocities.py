from matplotlib import pyplot as plt
import pandas as pd
import os

from matplotlib import pyplot as plt
import pandas as pd
import os
import numpy as np


plt.figure()
for source in ['mesa', 'kepler']:
    for M in ['15.0']:
        for E in ['1.0']:
            for Mix in ['5.0']:
                for Ni in [0.05]:
                    lum = pd.read_csv(os.path.join('data', 'example_velocities_snec',
                                                 'M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix)+'_R0_K0_'+source,'vel_photo.dat'),
                                    names=['time', 'veloc'], sep=r'\s+')
                    lum['time'] = lum['time'] / 86400
                    lum['veloc'] = lum['veloc'] / 100000000
                    plt.plot(lum['time'], lum['veloc'], label=source)
                    plt.legend()
                    plt.xlim(-5, 120)
                    plt.yticks(np.arange(0, 20, 5))
                    plt.ylim(0, 20)
                    plt.xlabel('time (days)')
                    plt.ylabel('Expansion velocity (10^3 km/s)', fontsize=11)
                    # plt.yscale('log')
                    # plt.xscale('log')
                    # plt.ylim(float(5*10**40), float(3.5*10**42))
                    plt.title('M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix))
                    plt.tight_layout()


plt.savefig(os.path.join('figures','velocities_M'+str(M)+'_Ni'+str(Ni)+'_E'+str(E)+'_Mix'+str(Mix)+'_R0_K0.png'))

plt.show()





plt.show()
