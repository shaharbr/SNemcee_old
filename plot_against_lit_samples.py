import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=30)     # fontsize of the axes title
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'




KSP_s50v = 0.3124
KSP_s50v_err = 0.006248
KSP_Vmax = -17.95999
KSP_Vmax_err = 0.39
KSP_v50 = -17.2099999
KSP_v50_err = 0.05
KSP_Ni_mass = 0.10
KSP_Ni_err = 0.01




SN2018hmx_param_Vband = pd.read_csv(r'results\SN2018hmx_param_results.csv')
SN2018hmx_Ni = pd.read_csv(r'results\Ni_results_SN2018hmx_BVgri.csv')

SN2018aad_param_Vband = pd.read_csv(r'results\SN2018aad_param_results.csv')
SN2018aad_Ni = pd.read_csv(r'results\Ni_results_SN2018aad_BVgri.csv')




def object_to_numeric(series):
    series = series.str.strip()
    series = pd.to_numeric(series)
    return series



valenti_s50v = pd.read_csv(r'data\valenti_result_s50_V2.tsv', sep='\t',
                           names=['Name', 's50V', 'e_s50V', 'NA1', 'NA2', 'NA3', 'NA4', 'NA5', ],
                           usecols=['Name', 's50V', 'e_s50V'])
valenti_s50v['e_s50V'] = object_to_numeric(valenti_s50v['e_s50V'])
valenti_s50v['s50V'] = 50 * valenti_s50v['s50V']
valenti_s50v['e_s50V'] = 50 * valenti_s50v['e_s50V']


valenti_bolometric = pd.read_csv(r'data\valenti2016_bolometric.tsv', sep='\t', header=51,
                                 usecols=['Name', 'MNi', 'e_MNi'])
valenti_bolometric.drop([0, 1], inplace=True)
for column in ['MNi', 'e_MNi']:
    valenti_bolometric[column] = object_to_numeric(valenti_bolometric[column])




valenti_vmag_50 = pd.read_csv(r'data\valenti_result_mag50_2.tsv', sep='\t',
                              names=['Name', 'vmag_50', 'e_vmag_50', 'jd', 'filter', 'NA1', 'NA2', 'NA3'],
                              usecols=['Name', 'vmag_50', 'e_vmag_50'])
valenti_vmag_50['e_vmag_50'] = object_to_numeric(valenti_vmag_50['e_vmag_50'])



valenti_vmax = pd.read_csv(r'data\valenti_result_mag_at_max.tsv', sep='\t',
                              names=['Name', 'Vmax', 'e_Vmax', 'NA1', 'NA2', 'NA3', 'NA4', 'NA5', 'NA6', 'NA7'],
                              usecols=['Name', 'Vmax', 'e_Vmax'])







for df in [valenti_s50v, valenti_vmax, valenti_vmag_50, valenti_bolometric]:
    df['Name'] = df['Name'].str.strip()
valenti_SN_param = valenti_s50v\
    .merge(valenti_vmax, on='Name', how='outer')\
    .merge(valenti_vmag_50, on='Name', how='outer')\
    .merge(valenti_bolometric, on='Name', how='outer')


fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(y=valenti_SN_param['Vmax'], yerr=valenti_SN_param['e_Vmax'],
            x=valenti_SN_param['s50V'], xerr=valenti_SN_param['e_s50V'],
            marker='.', fillstyle='none', Linestyle='None', color='k', label='Normal type II SNe', alpha=0.2)
ax.errorbar(y=SN2018hmx_param_Vband['Vmax'], yerr=SN2018hmx_param_Vband['e_Vmax'],
            x=SN2018hmx_param_Vband['s50V'], xerr=SN2018hmx_param_Vband['s50V_err'],
            marker='o', fillstyle='full', markersize=8, Linestyle='None', color='#FF3333', label='SN 2018hmx')

# ax.arrow(float(SN2018aad_param_Vband['s50V']), float(SN2018aad_param_Vband['Vmax']), 0.05, 0,
#          shape='full', color='#00C0C0', head_length=0.03, head_width=0.1)

ax.errorbar(y=SN2018aad_param_Vband['Vmax'], yerr=SN2018aad_param_Vband['e_Vmax'],
            x=SN2018aad_param_Vband['s50V'], xerr=SN2018aad_param_Vband['s50V_err'],
            marker='o', fillstyle='full', markersize=8, Linestyle='None', color='#00C0C0', label='SN 2018aad')
# ax.errorbar(y=KSP_Vmax, yerr=KSP_Vmax_err,
#             x=KSP_s50v, xerr=KSP_s50v_err,
#             marker='.', fillstyle='full', markersize=5, Linestyle='None', color='#CC6600', label='KSP-SN-2016kf',
#             alpha=0.6)
ax.invert_yaxis()
ax.set_title('Peak magnitude vs Decline rate\n')
ax.set_ylabel('Peak magnitude (V-band)')
ax.set_xlabel('Decline rate (V mag/50d)')
ax.set_xlim(left=-0.2)
ax.legend()
ax.tick_params(axis='both', which='major')
fig.savefig(r'figures\SN2018hmx_Vmax_s50V_against_valenti' + '.png')
fig.savefig(r'figures\SN2018hmx_Vmax_s50V_against_valenti' + '.svg')


fig, ax = plt.subplots(figsize=(9, 7))
ax.errorbar(y=valenti_SN_param['MNi'], yerr=valenti_SN_param['e_MNi'],
            x=valenti_SN_param['vmag_50'], xerr=valenti_SN_param['e_vmag_50'],
            marker='.', fillstyle='none', Linestyle='None', color='k', label='Normal type II SNe', alpha=0.2)


ax.errorbar(y=SN2018hmx_Ni['Ni_87A'], yerr=SN2018hmx_Ni['Ni_87A_sigma'],
            x=SN2018hmx_param_Vband['vmag_50'], xerr=SN2018hmx_param_Vband['vmag_50_err'],
            marker='o', fillstyle='full', markersize=7, Linestyle='None', color='#FF3333', capsize=2,
            label='SN 2018hmx\n by ratio to 1987A')

# ax.arrow(float(SN2018hmx_param_Vband['vmag_50']), float(SN2018hmx_Ni['Ni_mcmc']), 0, 0.05,
#          shape='full', color='#FF3333', head_length=0.03, head_width=0.1)
ax.errorbar(y=SN2018hmx_Ni['Ni_mcmc'], yerr=[SN2018hmx_Ni['Ni_mcmc_16perc'], SN2018hmx_Ni['Ni_mcmc_84perc']],
            x=SN2018hmx_param_Vband['vmag_50'], xerr=SN2018hmx_param_Vband['vmag_50_err'],
            marker='o', fillstyle='full', markersize=8, Linestyle='None', color='#911310', capsize=2,
            label='SN 2018hmx\n by MCMC')

# ax.arrow(float(SN2018hmx_param_Vband['vmag_50']), 0.170639494, 0, 0.05,
#          shape='full', color='#911310', head_length=0.03, head_width=0.1)
ax.errorbar(y=SN2018aad_Ni['Ni_87A'], yerr=SN2018aad_Ni['Ni_87A_sigma'],
            x=SN2018aad_param_Vband['vmag_50'], xerr=SN2018aad_param_Vband['vmag_50_err'],
            marker='o', fillstyle='full', markersize=7, Linestyle='None', color='#2db5a1', capsize=2,
            label='SN 2018aad\n by ratio to 1987A')


ax.errorbar(y=SN2018aad_Ni['Ni_mcmc'], yerr=[SN2018aad_Ni['Ni_mcmc_16perc'], SN2018aad_Ni['Ni_mcmc_84perc']],
            x=SN2018aad_param_Vband['vmag_50'], xerr=SN2018aad_param_Vband['vmag_50_err'],
            marker='o', fillstyle='full', markersize=8, Linestyle='None', color='#166156', capsize=2,
            label='SN 2018aad\n by MCMC')


# ax.errorbar(y=KSP_Ni_mass, yerr=KSP_Ni_err,
#             x=KSP_v50, xerr=KSP_v50_err,
#             marker='.', fillstyle='full', markersize=5, Linestyle='None', color='#CC6600', label='KSP-SN-2016kf', alpha=0.6)

ax.invert_xaxis()
ax.set_yscale('log')
ax.set_title('Ni mass vs day 50 magnitude\n')
ax.set_ylabel('Ni mass ($M_{\odot}$)')
ax.set_xlabel('V band absolute mag at day 50')
ax.set_xlim(left=-13, right=-19)
ax.set_ylim(top=0.3)
ax.legend(loc='upper left')
ax.tick_params(axis='y', which='major', labelsize=24)
fig.savefig(r'figures\SN2018hmx_vmag50_Ni_against_valenti' + '.png')
fig.savefig(r'figures\SN2018hmx_vmag50_Ni_against_valenti' + '.svg')


gutierrez_veloc = pd.read_csv(r'data\Gutierrez_2017_apjaa8f52t8_ascii.txt', sep='\t', header=2,
                           usecols=['Epoch', r'${{\rm{H}}}_{\alpha }$.1', r'${{\rm{H}}}_{\beta }$',
                                    'Fe II lambda5169'])
gutierrez_veloc.rename(columns={r'${{\rm{H}}}_{\alpha }$.1': 'Halpha',
                             r'${{\rm{H}}}_{\beta }$': 'Hbeta', 'Fe II lambda5169': 'FeII 5169'}, inplace=True)
gutierrez_veloc.drop([0, 24], inplace=True)

gutierrez_veloc['Halpha'], gutierrez_veloc['e_Halpha'] = gutierrez_veloc['Halpha'].str.split('or', 2).str
gutierrez_veloc['Hbeta'], gutierrez_veloc['e_Hbeta'] = gutierrez_veloc['Hbeta'].str.split('or', 2).str
gutierrez_veloc['FeII 5169'], gutierrez_veloc['e_FeII 5169'] = gutierrez_veloc['FeII 5169'].str.split('or', 2).str



for column in gutierrez_veloc.keys():
    gutierrez_veloc[column].replace(regex=True, inplace=True, to_replace=r'\+|\-|\s|cdots', value='')
    gutierrez_veloc[column] = pd.to_numeric(gutierrez_veloc[column])



colors = {'Halpha': '#1b9e77', 'Hbeta': '#7570b3', 'FeII 5169': '#d95f02'}

SN2018hmx_veloc = pd.read_csv(r'results\sN2018hmx_expansion_velocities.csv')
iPTF14hls_veloc = pd.read_csv(r'results\SNiPTF14hls_expansion_velocities.csv')
SN2018aad_veloc = pd.read_csv(r'results\SN2018aad_expansion_velocities.csv')
# TODO correct 14hls time from discovery

SN_veloc_dict = {'SN 2018hmx': SN2018hmx_veloc, 'iPTF14hls': iPTF14hls_veloc, 'SN 2018aad': SN2018aad_veloc}
formatting_dict = {'SN 2018hmx': ['o', 'full', 'None'],
                   'iPTF14hls': ['s', 'full', 'None'],
                   'SN 2018aad': ['^', 'full', 'None']}




def plot_velocities(SN_veloc_dict, gutierrez_veloc, line_list, formatting_dict):
    line_names_formatted = {'Halpha': r'H$\alpha$', 'Hbeta': r'H$\beta$', 'FeII 5169': r'FeII 5169$\AA$'}
    line_num = len(line_list)
    fig2, axes = plt.subplots(line_num, 1, sharex='col', figsize=(10, 10))
    for i in range(len(line_list)):
        line_name = line_list[i]
        # gutierrez sample range
        axes[i].fill_between(gutierrez_veloc['Epoch'],
                        gutierrez_veloc[line_name] - gutierrez_veloc['e_' + line_name],
                        gutierrez_veloc[line_name] + gutierrez_veloc['e_' + line_name],
                        linestyle='--',
                        color='gray',
                        alpha=0.1)
        # SN-of-interest mean + error bars
        max_y = []
        for SN_name in SN_veloc_dict.keys():
            if SN_name == 'SN 2018hmx':
                color = '#1b9e77'
            elif SN_name == 'SN 2018aad':
                color = '#d95f02'
            elif SN_name == 'iPTF14hls':
                color = '#7570b3'
            #     color = colors[line_name]
            # else:
            #     color = 'gray'

            SN_veloc = SN_veloc_dict[SN_name]
            linedf = SN_veloc.loc[SN_veloc['line'] == line_name]
            axes[i].errorbar(x=linedf['t_from_discovery'], y=linedf['absorption_mean_velocity'],
                        yerr=linedf['absorption_std_velocity'],
                        label=SN_name,
                        marker=formatting_dict[SN_name][0],
                        fillstyle=formatting_dict[SN_name][1],
                        linestyle=formatting_dict[SN_name][2],
                        color=color)
            max_y = np.append(max_y, np.max(linedf['absorption_mean_velocity']))
        max_y = np.max(max_y)
        ticks = np.arange(0, max_y + 3000, 2500)
        axes[i].set_yticks(ticks)
        axes[i].set_ylim(top=np.max(ticks)*1.15)
        axes[i].text(160, np.max(ticks)*0.95, line_names_formatted[line_name], horizontalalignment='center')
        # gutierrez sample mean
        axes[i].plot(gutierrez_veloc['Epoch'], gutierrez_veloc[line_name],
                     # label='Gutierrez at al. mean $\pm$std',
                     label='Normal type II SNe',
                     marker='None',
                     linestyle='--',
                     color='gray',
                     alpha=0.4)
        # axes[i].set_title(line_names_formatted[line_name])


    # axes[0].errorbar(x=93, y=7960, yerr=60,
    #                  label='KSP-SN-2016kf',
    #                  marker='*',
    #                  markersize=10,
    #                  fillstyle='none',
    #                  linestyle='None',
    #                  color='gray')

    # axes[1].errorbar(x=93, y=5860, yerr=45,
    #                  label='KSP-SN-2016kf',
    #                  marker='*',
    #                  markersize=10,
    #                  fillstyle='none',
    #                  linestyle='None',
    #                  color='gray')
    #
    # axes[2].errorbar(x=93, y=3192, yerr=63,
    #                  label='KSP-SN-2016kf',
    #                  marker='*',
    #                  markersize=10,
    #                  fillstyle='none',
    #                  linestyle='None',
    #                  color='gray')


    for ax in axes:
        ax.legend(loc='upper right')
        ax.set_ylim(bottom=0)
        ax.set_xlim(-5, 350)
        ax.tick_params(axis='both', which='major')
    axes[1].set_ylabel('Expansion velocity (km/s)')
    plt.xlabel('Days from discovery')
    # plt.suptitle('Expansion velocity over time')
    fig2.savefig(r'figures\SN2018hmx_SN2018aad_iPTF14hls_absorption_velocities_against_gutierrez' + '.png')
    fig2.savefig(r'figures\SN2018hmx_SN2018aad_iPTF14hls_absorption_velocities_against_gutierrez' + '.svg')


plot_velocities(SN_veloc_dict, gutierrez_veloc, ['FeII 5169', 'Halpha'], formatting_dict)

plt.show()