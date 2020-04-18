import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=30)     # fontsize of the axes title
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



velocities_day_50 = pd.read_csv(r'data\Gutierrez_Velocities_day_50.csv')
A14_vmag = pd.read_csv(r'data\A14_vmag_end.csv')
A14_vmag = A14_vmag.merge(velocities_day_50)


SN2018hmx_param_Vband = pd.read_csv(r'results\SN2018hmx_param_results.csv')
SN2018aad_param_Vband = pd.read_csv(r'results\SN2018aad_param_results.csv')



def find_closest_to_50(df, minORmax):
    daylist = df['t_from_discovery']
    if minORmax == 'min':
        daylist = daylist[daylist > 50]
        minday = np.min(daylist)
    elif minORmax == 'max':
        daylist = daylist[daylist < 50]
        minday = np.max(daylist)
    return minday

def get_v_closest_to_50(df, minORmax):
    day = find_closest_to_50(df, minORmax)
    df = df.loc[(df['line'] == 'FeII 5169') & (df['t_from_discovery'] == day)]
    v = df['absorption_mean_velocity']
    d_v = df['absorption_std_velocity']
    return v, d_v, day

SN2018hmx_veloc = pd.read_csv(r'results\sN2018hmx_expansion_velocities.csv')
SN18hmx_v, SN18hmx_dv, SN18hmx_day = get_v_closest_to_50(SN2018hmx_veloc, 'max')

iPTF14hls_veloc = pd.read_csv(r'results\SNiPTF14hls_expansion_velocities.csv')
iPTF14hls_v, iPTF14hls_dv, iPTF14hls_day = get_v_closest_to_50(iPTF14hls_veloc, 'max')

SN2018aad_veloc = pd.read_csv(r'results\SN2018aad_expansion_velocities.csv')
SN18aad_v, SN18aad_dv, SN18aad_day = get_v_closest_to_50(SN2018aad_veloc, 'min')


fig, ax = plt.subplots(figsize=(8, 8))
ax.errorbar(y=A14_vmag['M_end'], yerr=A14_vmag['d_M_end'],
            x=A14_vmag['Fe_v_day50'], xerr=A14_vmag['d_Fe_v_day50'],
            # marker='o', fillstyle='none', Linestyle='None', color='k', label='A14 sample M_end\n(30 days before tPT)', alpha=0.2)
            marker='o', fillstyle='none', Linestyle='None', color='k', label='Normal type II SNe', alpha=0.2)



# ax.arrow(float(SN18hmx_v), float(SN2018hmx_param_Vband['vmag_50']), 100, 0,
#          shape='full', color='#FF3333', head_length=50, head_width=0.1)

ax.errorbar(y=SN2018hmx_param_Vband['vmag_50'], yerr=SN2018hmx_param_Vband['vmag_50_err'],
            x=SN18hmx_v, xerr=SN18hmx_dv,
            marker='o', markersize=8, fillstyle='full', Linestyle='None', color='#FF3333', label='SN18hmx day '+str(int(SN18hmx_day)), alpha=0.7)

ax.arrow(float(SN18aad_v), float(SN2018aad_param_Vband['vmag_50']), 300, 0,
         shape='full', color='#00C0C0', head_length=70, head_width=0.12)

ax.errorbar(y=SN2018aad_param_Vband['vmag_50'], yerr=SN2018aad_param_Vband['vmag_50_err'],
            x=SN18aad_v, xerr=SN18aad_dv,
            marker='o', markersize=8, fillstyle='full', Linestyle='None', color='#00C0C0', label='SN18aad day '+str(int(SN18aad_day)), alpha=0.7)

ax.errorbar(y=-19, yerr=SN2018aad_param_Vband['vmag_50_err'],
            x=iPTF14hls_v, xerr=iPTF14hls_dv,
            marker='o', markersize=8, fillstyle='full', Linestyle='None', color='#00C0C0', label='SN18aad day '+str(int(SN18aad_day)), alpha=0.7)


ax.legend()
ax.invert_yaxis()
ax.set_ylim(top=-20.5)
ax.set_xlim(right=6500)
ax.set_title('V-band magnitude vs. \n expansion velocity around day 50\n')
# ax.set_xlabel('v(Fe 5169) around day 50 (km/s)')
# ax.set_ylabel('Absolute V magnitude at day 50 (mag)')

ax.set_xlabel('Expansion velocity (day 50, km/s)\n', fontsize=24)
ax.set_ylabel('\nMagnitude (day 50 V-band mag)', fontsize=24)

plt.show()