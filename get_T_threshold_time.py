import numpy as np
import pandas as pd
import os
import corner
from matplotlib import pyplot as plt


plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=20)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



T_thresh = 10 ** 3.75

Mzams_range = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
Ni_range = [0.02, 0.07, 0.12]
E_final_range = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]
Mix_range = [2.0, 5.0, 8.0]

M_list = []
Ni_list = []
E_list = []
Mix_list = []
T_list = []


for M in Mzams_range:
    for Ni in Ni_range:
        for E in E_final_range:
            for Mix in Mix_range:
                M_list += [M]
                Ni_list += [Ni]
                E_list += [E]
                Mix_list += [Mix]
                model = 'M' + str(M) + '_Ni' + str(Ni) + '_E' + str(E) + '_Mix' + str(Mix) + '_R0_K0'
                temps = pd.read_csv(os.path.join('all_temp_rad_data', model, 'T_eff.dat'), names=['time', 'temp'], sep=r'\s+')
                max_temp_below_Tthresh = np.max(temps['temp'].loc[temps['temp'] < T_thresh])
                time_thresh = float(temps['time'].loc[temps['temp'] == max_temp_below_Tthresh]) / 86400
                T_list += [time_thresh]


all_array = np.array([M_list, Ni_list, E_list, Mix_list, T_list]).T

print(all_array)
labels = ['Mzams', 'Ni', 'E', 'Mix', 'T_threshold']
corner_range = [1., 1., 1., 1., 1.]
f_corner = corner.corner(all_array, labels=labels, range=corner_range)
plt.tight_layout()

f_corner.savefig(os.path.join('figures', 'corner_plot_T_threshold.png'))





