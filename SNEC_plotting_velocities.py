import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import os

data_dir = os.path.join('data', 'example_velocities_snec')

for Mzams in ['9.0', '15.0', '17.0']:
    name = 'M'+Mzams+'_Ni0.02_E0.7_Mix2.0_R600_K0.001'
    snec_model = pd.read_csv(os.path.join(data_dir, name, 'vel_photo.dat'),
                                     names=['t_from_discovery', 'veloc'], sep=r'\s+')
    # convert sec to days
    snec_model['t_from_discovery'] = snec_model['t_from_discovery'] / 86400
    # convery cm/s to km/s
    snec_model['veloc'] = snec_model['veloc'] / 100000
    time_col = snec_model['t_from_discovery']
    snec_model = snec_model['veloc']
    snec_model = np.interp(np.arange(0, 400, 10), time_col, snec_model)
    plt.plot(np.arange(0, 400, 10), snec_model, label=Mzams)
    plt.legend()



data_dir = os.path.join('data', 'example_velocities_snec')

for Mzams in ['9.0', '15.0']:
    name = 'M'+Mzams+'_Ni0.02_E0.7_Mix2.0_R600_K0.001'
    snec_model = pd.read_csv(os.path.join(data_dir, name, 'lum_observed.dat'),
                             names=['t_from_discovery', 'Lum'], sep=r'\s+')
    time_col = snec_model['t_from_discovery'] / 86400
    snec_model = snec_model['Lum']
    snec_model = np.interp(np.arange(0, 400, 1), time_col, snec_model)
    plt.plot(np.arange(0, 400, 1), snec_model, label=Mzams)
    plt.legend()



plt.show()
