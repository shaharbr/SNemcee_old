import pandas as pd
from matplotlib import pyplot as plt

SN2018hmx_param_Vband = pd.read_csv(r'results\SN2018hmx_param_results.csv')
SN2018hmx_Ni = pd.read_csv(r'results\SN2018hmx_Ni_results.csv')


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





# TODO Ni mass against Vmag at day 50
# TODO error bars for Vmax too small?
# TODO how to calculate std for SV50 and other params?
# TODO use mag at 50 from valenti
# TODO check valenti method for bolometric calculation
# TODO add ATLAS data


for df in [valenti_s50v, valenti_vmax, valenti_vmag_50, valenti_bolometric]:
    df['Name'] = df['Name'].str.strip()
valenti_SN_param = valenti_s50v\
    .merge(valenti_vmax, on='Name', how='outer')\
    .merge(valenti_vmag_50, on='Name', how='outer')\
    .merge(valenti_bolometric, on='Name', how='outer')


fig, ax = plt.subplots()
ax.errorbar(y=valenti_SN_param['Vmax'], yerr=valenti_SN_param['e_Vmax'],
            x=valenti_SN_param['s50V'], xerr=valenti_SN_param['e_s50V'],
            marker='.', fillstyle='none', Linestyle='None', color='k', label='Valenti SNe')
ax.errorbar(y=SN2018hmx_param_Vband['Vmax'], yerr=SN2018hmx_param_Vband['e_Vmax'],
            x=SN2018hmx_param_Vband['s50V'], xerr=SN2018hmx_param_Vband['s50V_err'],
            marker='.', fillstyle='none', Linestyle='None', color='orange', label='SN2018hmx')
# TODO need errorbars for 18hmx
ax.invert_yaxis()
ax.set_title('Vmax vs s50V', fontsize=16)
ax.set_ylabel('Vmax', size=12)
ax.set_xlabel('s50V', size=12)
ax.set_xlim(left=-0.2)
ax.legend()
ax.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(r'figures\SN2018hmx_Vmax_s50V_against_valenti' + '.png')






fig, ax = plt.subplots()
ax.errorbar(y=valenti_SN_param['MNi'], yerr=valenti_SN_param['e_MNi'],
            x=valenti_SN_param['vmag_50'], xerr=valenti_SN_param['e_vmag_50'],
            marker='.', fillstyle='none', Linestyle='None', color='k', label='Valenti SNe')
# TODO get the 18hmx values from the table

ax.errorbar(y=SN2018hmx_Ni['Ni_mass'], yerr=[SN2018hmx_Ni['Ni_lower'], SN2018hmx_Ni['Ni_upper']],
            x=SN2018hmx_param_Vband['vmag_50'], xerr=SN2018hmx_param_Vband['vmag_50_err'],
            marker='.', fillstyle='none', Linestyle='None', color='orange', label='SN2018hmx')
# TODO need errorbars for 18hmx
ax.invert_xaxis()
ax.set_yscale('log')
ax.set_title('Ni mass vs V band absolute mag at day 50', fontsize=16)
ax.set_ylabel('Ni mass (solar masses)', size=12)
ax.set_xlabel('V band absolute mag at day 50', size=12)
ax.legend()
ax.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(r'figures\SN2018hmx_vmag50_Ni_against_valenti' + '.png')





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
SNiPTF14hls_veloc = pd.read_csv(r'results\SNiPTF14hls_expansion_velocities.csv')
SN2018aad_veloc = pd.read_csv(r'results\SN2018aad_expansion_velocities.csv')
# TODO correct 14hls time from discovery


fig, ax = plt.subplots(figsize=(12, 6))

for line_name in colors.keys():
    ax.fill_between(gutierrez_veloc['Epoch'],
                    gutierrez_veloc[line_name] - gutierrez_veloc['e_' + line_name],
                    gutierrez_veloc[line_name] + gutierrez_veloc['e_' + line_name],
                    linestyle='--',
                    color=colors[line_name],
                    alpha=0.2)

    ax.plot(gutierrez_veloc['Epoch'], gutierrez_veloc[line_name],
            label='Gutierrez ' + line_name + r' mean $\pm$std',
            marker='None',
            linestyle='--',
            color=colors[line_name],
            alpha=0.6)

    linedf = SN2018hmx_veloc.loc[SN2018hmx_veloc['line'] == line_name]
    ax.errorbar(x=linedf['t_from_discovery'], y=linedf['absorption_mean_velocity'],
                yerr=linedf['absorption_std_velocity'],
                label='SN2018hmx ' + line_name,
                marker='s',
                fillstyle='full',
                linestyle='None',
                color=colors[line_name])

    linedf = SNiPTF14hls_veloc.loc[SNiPTF14hls_veloc['line'] == line_name]
    ax.errorbar(x=linedf['t_from_discovery'], y=linedf['absorption_mean_velocity'],
                yerr=linedf['absorption_std_velocity'],
                label='SNiPTF14hls ' + line_name,
                marker='s',
                fillstyle='none',
                linestyle='None',
                color=colors[line_name])

    linedf = SN2018aad_veloc.loc[SN2018aad_veloc['line'] == line_name]
    ax.errorbar(x=linedf['t_from_discovery'], y=linedf['absorption_mean_velocity'],
                yerr=linedf['absorption_std_velocity'],
                label='SN2018aad ' + line_name,
                marker='^',
                fillstyle='none',
                linestyle='None',
                color=colors[line_name])



ax.set_title('Expansion velocity over time', fontsize=16)
ax.set_ylabel('Expansion velocity (km/s)', size=12)
ax.set_xlabel('Days from discovery', size=12)
ax.legend(ncol=3)
ax.set_ylim(bottom=0)
ax.tick_params(axis='both', which='major', labelsize=14)
fig.savefig(r'figures\SN2018hmx_SN2018aad_SNiPTF14hls_absorption_velocities_against_gutierrez' + '.png')



plt.show()