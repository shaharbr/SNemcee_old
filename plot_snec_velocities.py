from matplotlib import pyplot as plt
import pandas as pd
import os

snec_model = pd.read_csv(os.path.join('..', 'all_veloc_data',
                                 'M11.0_Ni0.02_E1.0_Mix2.0_R600_K0.001','vel_photo.dat'),
                    names=['t_from_discovery', 'veloc'], sep=r'\s+')
# convert sec to days
snec_model['t_from_discovery'] = snec_model['t_from_discovery'] / 86400
# convery cm/s to km/s
snec_model['veloc'] = snec_model['veloc'] / 100000
time_col = snec_model['t_from_discovery']
snec_model = snec_model['veloc']

plt.figure()
plt.plot(time_col, snec_model)
plt.title('Expansion velocity over time\n'
          'SNEC model: M11.0_Ni0.02_E1.0_Mix2.0_R600_K0.001', fontsize=12)
plt.ylabel('Expansion velocity (km/s)', fontsize=11)
plt.xlabel('Time (days)', fontsize=10)
plt.savefig(os.path.join('results', 'SNEC_M11.0_Ni0.02_E1.0_Mix2.0_R600_K0.001_velocities_400d.png'))



snec_model = pd.read_csv(os.path.join('..', 'all_lum_data',
                                 'M11.0_Ni0.02_E1.0_Mix2.0_R600_K0.001','lum_observed.dat'),
                    names=['t_from_discovery', 'lum'], sep=r'\s+')
# convert sec to days
snec_model['t_from_discovery'] = snec_model['t_from_discovery'] / 86400
# convery cm/s to km/s
snec_model['lum'] = snec_model['lum']
time_col = snec_model['t_from_discovery']
snec_model = snec_model['lum']

plt.figure()
plt.plot(time_col, snec_model)
plt.title('Bolometric luminosity over time\n'
          'SNEC model: M11.0_Ni0.02_E1.0_Mix2.0_R600_K0.001', fontsize=12)
plt.ylabel('Bolometric luminosity (erg/s)', fontsize=11)
plt.yscale('log')
plt.ylim(8.4e+39, 5.2e+42)
plt.xlabel('Time (days)', fontsize=10)
plt.savefig(os.path.join('results', 'SNEC_M11.0_Ni0.02_E1.0_Mix2.0_R600_K0.001_lum_400d.png'))



plt.show()
