import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from os import path

# plt.rc('font', size=20)          # controls default text sizes
# plt.rc('axes', titlesize=20)     # fontsize of the axes title
# plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=8)    # legend fontsize
# plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'


# def calc_chi_square(data, model):
#     data_after_20 = data.loc[data['t_from_discovery'] > 20]
#     interp_model_lum = np.interp(data_after_20['t_from_discovery'],
#                                     model['t_from_discovery'],
#                                     model['Lum'])
#     chisq = np.sum(((data_after_20['Lum'] - interp_model_lum) /
#                     np.mean([data_after_20['dLum0'], data_after_20['dLum1']])) ** 2)
#     chisq_reduced = chisq / (len(interp_model_lum) - 1)
#     return chisq_reduced


def calc_chi_square_sampled(data, model):
    sampling = np.arange(20, 170, 5)
    data_sampled = np.interp(sampling, data['t_from_discovery'], data['Lum'])
    data_err_sampled = np.interp(sampling, data['t_from_discovery'], data['dLum0'])
    model_sampled = np.interp(sampling, model['t_from_discovery'], model['Lum'])
    chisq = np.sum(((data_sampled - model_sampled) /
                    data_err_sampled) ** 2)
    chisq_reduced = chisq / (len(sampling) - 1)
    # plt.figure()
    # plt.plot(sampling, data_sampled, marker='o')
    # plt.plot(sampling, model_sampled, marker='o')
    # plt.yscale('log')
    return chisq_reduced



# obs_mag_016 = pd.read_csv(r'Data_Ni012_E120\magnitudes.dat',
#                           names=['t_from_discovery', 'unknown', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'],
#                           sep=r'\s+')
# print(obs_mag_016)
# obs_mag_016['t_from_discovery'] = obs_mag_016['t_from_discovery'] /1440



blackbody_data = pd.read_csv(r'results\blackbody_results_18hmx.csv')

# convert watt to erg/s
blackbody_data['Lum'] = blackbody_data['Lum'] * 10**7
blackbody_data['dLum0'] = blackbody_data['dLum0'] * 10**7
blackbody_data['dLum1'] = blackbody_data['dLum1'] * 10**7


sims =['M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)\
                   + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
                    for Mzams in [13.0, 19.0]
                    for Ni_mass in [0.10, 0.16]
                    for E_final in [1.5, 2.4]
                    for Ni_boundary in [3.0]
                    for R_CSM in [600, 2400]
                    for K_CSM in [0.001, 90]]
print(sims)


def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]





lums = {}



# colors = [(r, g, b) for r,g,b in zip([1., .7, .8, .9, .8, .7, .5, .1, .0, .1, .1, .1, .1, .1, .1, .1, .0, .1, .5, .7, .8],
#                                      [.0, .2, .5, .7, .8, .9, .8, .9, 1., .7, .8, .9, .8, .7, .5, .1, .0, .2, .1, .1, .1],
#                                      [.0, .1, .1, .1, .1, .1, .1, .1, .0, .2, .5, .7, .8, .9, .8, .9, 1., .7, .8, .9, .7],
#                                      )]

# colors = [(r, g, b) for r,g,b in zip([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., .9, .8, .7, .6, .5, .4, .3, .2, .1, .0],
#                                      [.1, .1, .1, .1, .1, .1, .1, .5, .5, .5, .5, .5, .5, .5, .9, .9, .9, .9, .9, .9, .9],
#                                      [.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
#                                      )]





titles = sims


cols = 4
rows = 8

num_graphs = cols * rows
third = int(num_graphs/3) + 1

colors = [(r, g, b) for r,g,b in zip(list(np.arange(1.0, 0.0, -1/third)) + list(np.arange(0.0, 1.0, .8/third))+ [0.0] * third,
                                     [0.0] * third + list(np.arange(0.0, 1.0, .8/third)) + [0.0] * third + list(np.arange(1.0, 0.0, -0.8/third)),
                                     list(np.arange(0.0, 1.0, .8/third)) + list(np.arange(1.0, 0.0, -0.8/third)) + [0.0] * third
                                     )]

print(colors)


fig1, axs = plt.subplots(rows, cols, figsize=(10, 20), constrained_layout=True, sharey=True, sharex=True)
fig2, ax2 = plt.subplots(figsize=(10,8), constrained_layout=True)

plt.ylabel('bolometric luminosity')
plt.xlabel('days')
plt.title('wind parameters')


axs = trim_axs(axs, len(sims))
for ax1, sim, color, title in zip(axs, sims, colors, titles):
        lums[sim] = pd.read_csv(path.join('..\Documents', 'all_lum_data', sim, 'lum_observed.dat'), names=['t_from_discovery', 'Lum'], sep=r'\s+')
        lums[sim]['t_from_discovery'] = lums[sim]['t_from_discovery'] /86400

        chisq = calc_chi_square_sampled(blackbody_data, lums[sim])

        for ax in [ax1, ax2]:
            ax.plot(lums[sim]['t_from_discovery'],
                    lums[sim]['Lum'],
                    marker='.', markersize=2, Linestyle='None', label=title, alpha=1, color=color)

        ax1.errorbar(x=blackbody_data['t_from_discovery'],
                    y=blackbody_data['Lum'],
                    yerr=[blackbody_data['dLum0'], blackbody_data['dLum1']],
                    label='SN 2018hmx', marker='.',  markersize=5, fillstyle='full', linestyle='None', color='k')

        ax1.set_yscale('log')
        ax1.set_ylim(2.8e+41, 1.3e+43)
        ax1.set_title('chi_sq_red = '+str(int(chisq)), fontsize=14)
        # ax1.set_ylabel('bolometric luminosity')
        # ax1.set_xlabel('days')
        # plt.title('wind parameters')
        ax1.legend()




ax2.errorbar(x=blackbody_data['t_from_discovery'],
            y=blackbody_data['Lum'],
            yerr=[blackbody_data['dLum0'], blackbody_data['dLum1']],
            label='SN 2018hmx', marker='.',  markersize=5, fillstyle='full',linestyle='None', color='k')
ax2.set_yscale('log')
ax2.set_ylim(2.8e+41, 1.3e+43)
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





