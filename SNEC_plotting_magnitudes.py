import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

# plt.rc('font', size=20)          # controls default text sizes
# plt.rc('axes', titlesize=20)     # fontsize of the axes title
# plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=8)    # legend fontsize
# plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



# obs_mag_016 = pd.read_csv(r'Data_Ni012_E120\magnitudes.dat',
#                           names=['t_from_discovery', 'unknown', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'],
#                           sep=r'\s+')
# print(obs_mag_016)
# obs_mag_016['t_from_discovery'] = obs_mag_016['t_from_discovery'] /1440



blackbody_data = pd.read_csv(r'results\blackbody_results.csv')

# convert watt to erg/s
blackbody_data['Lum'] = blackbody_data['Lum'] * 10**7
blackbody_data['dLum0'] = blackbody_data['dLum0'] * 10**7
blackbody_data['dLum1'] = blackbody_data['dLum1'] * 10**7

sims = [
        'Ni014_E240_Mix1_Mzams15',
        'Ni014_E240_Mix1_Mzams18',
        'Ni014_E240_Mix1_Mzams21',

        'Ni014_E240_Mix3_Mzams15',
        'Ni014_E240_Mix3_Mzams18',
        'Ni014_E240_Mix3_Mzams21',

        'Ni014_E270_Mix1_Mzams15',
        'Ni014_E270_Mix1_Mzams18',
        'Ni014_E270_Mix1_Mzams21',

        'Ni014_E270_Mix3_Mzams15',
        'Ni014_E270_Mix3_Mzams18',
        'Ni014_E270_Mix3_Mzams21',
        ]


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]


lums = {}

colors = [(r, g, b) for r,g,b in zip([1., .7, .8, .9, .8, .7, .5, .1, .0, .1, .1, .1, .1, .1, .1, .1, .0, .1, .5, .7, .8],
                                     [.0, .2, .5, .7, .8, .9, .8, .9, 1., .7, .8, .9, .8, .7, .5, .1, .0, .2, .1, .1, .1],
                                     [.0, .1, .1, .1, .1, .1, .1, .1, .0, .2, .5, .7, .8, .9, .8, .9, 1., .7, .8, .9, .7],
                                     )]

# colors = [(r, g, b) for r,g,b in zip([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., .9, .8, .7, .6, .5, .4, .3, .2, .1, .0],
#                                      [.1, .1, .1, .1, .1, .1, .1, .5, .5, .5, .5, .5, .5, .5, .9, .9, .9, .9, .9, .9, .9],
#                                      [.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
#                                      )]





titles = [
        'Ni=0.14 E=2.4 Ni_mix=1 Mzams=15',
        'Ni=0.14 E=2.4 Ni_mix=1 Mzams=18',
        'Ni=0.14 E=2.4 Ni_mix=1 Mzams=21',

        'Ni=0.14 E=2.4 Ni_mix=3 Mzams=15',
        'Ni=0.14 E=2.4 Ni_mix=3 Mzams=18',
        'Ni=0.14 E=2.4 Ni_mix=3 Mzams=21',

        'Ni=0.14 E=2.7 Ni_mix=1 Mzams=15',
        'Ni=0.14 E=2.7 Ni_mix=1 Mzams=18',
        'Ni=0.14 E=2.7 Ni_mix=1 Mzams=21',

        'Ni=0.14 E=2.7 Ni_mix=3 Mzams=15',
        'Ni=0.14 E=2.7 Ni_mix=3 Mzams=18',
        'Ni=0.14 E=2.7 Ni_mix=3 Mzams=21',
        ]


cols = 3
rows = 4

fig1, axs = plt.subplots(rows, cols, figsize=(10, 20), constrained_layout=True, sharey=True, sharex=True)
fig2, ax2 = plt.subplots(figsize=(10,8), constrained_layout=True)



axs = trim_axs(axs, len(sims))
for ax1, sim, color, title in zip(axs, sims, colors, titles):
        lums[sim] = pd.read_csv(sim+r'\lum_observed.dat', names=['t_from_discovery', 'Lum'], sep=r'\s+')
        lums[sim]['t_from_discovery'] = lums[sim]['t_from_discovery'] /86400
        for ax in [ax1, ax2]:
            ax.plot(lums[sim]['t_from_discovery'],
                    lums[sim]['Lum'],
                    marker='.', markersize=2, Linestyle='None', label=title, alpha=1, color=color)

        ax1.errorbar(x=blackbody_data['t_from_discovery'],
                    y=blackbody_data['Lum'],
                    yerr=[blackbody_data['dLum0'], blackbody_data['dLum1']],
                    label='SN 2018hmx', marker='.',  markersize=5, fillstyle='full', linestyle='None', color='k')

        ax1.set_yscale('log')
        # ax1.set_ylabel('bolometric luminosity')
        # ax1.set_xlabel('days')
        # ax1.set_title(title, fontsize=14)
        ax1.legend()


ax2.errorbar(x=blackbody_data['t_from_discovery'],
            y=blackbody_data['Lum'],
            yerr=[blackbody_data['dLum0'], blackbody_data['dLum1']],
            label='SN 2018hmx', marker='.',  markersize=5, fillstyle='full',linestyle='None', color='k')
ax2.set_yscale('log')
ax2.set_ylabel('bolometric luminosity')
ax2.set_xlabel('time from discovery (days)')
ax2.set_title('Explosion energy - overlay all', fontsize=14)
ax2.legend()


plt.show()


# lum_observed_Ni012_E120 = pd.read_csv(r'Data_Ni012_E120\lum_observed.dat', names=['t_from_discovery', 'Lum'],sep=r'\s+')
# lum_observed016 = pd.read_csv(r'SNEC_0.16NiMass\lum_observed.dat', names=['t_from_discovery', 'Lum'],sep=r'\s+')
# lum_observed016_2excis = pd.read_csv(r'SNEC_0.16NiMass_2Mexcised\lum_observed.dat', names=['t_from_discovery', 'Lum'],sep=r'\s+')
# lum_observed005 = pd.read_csv(r'SNEC_deafult_15MsRSG\lum_observed.dat', names=['t_from_discovery', 'Lum'],sep=r'\s+')
# convert seconds to days
# lum_observed016['t_from_discovery'] = lum_observed016['t_from_discovery'] /86400
# lum_observed005['t_from_discovery'] = lum_observed005['t_from_discovery'] /86400
# lum_observed016_2excis['t_from_discovery'] = lum_observed016_2excis['t_from_discovery'] /86400
# lum_observed_Ni012_E120['t_from_discovery'] = lum_observed_Ni012_E120['t_from_discovery'] /86400





# ax.plot(lum_observed_Ni012_E120['t_from_discovery'],
#         lum_observed_Ni012_E120['Lum'],
#         marker='.', fillstyle='none', Linestyle='None', color='red', label='0.16 Ms', alpha=0.4)
#
# ax.plot(lum_observed016['t_from_discovery'],
#         lum_observed016['Lum'],
#         marker='.', fillstyle='none', Linestyle='None', color='orange', label='0.16 Ms', alpha=0.4)
#
# ax.plot(lum_observed016_2excis['t_from_discovery'],
#         lum_observed016_2excis['Lum'],
#         marker='.', fillstyle='none', Linestyle='None', color='purple', label='0.16 Ms 2 excised', alpha=0.2)
#
# ax.plot(lum_observed005['t_from_discovery'],
#         lum_observed005['Lum'],
#         marker='.', fillstyle='none', Linestyle='None', color='green', label='0.05 Ms', alpha=0.4)

# for filter in  ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
#     ax.plot(obs_mag_016['t_from_discovery'],
#             obs_mag_016[filter],
#             marker='.', fillstyle='none', Linestyle='None', label=filter, alpha=1)





