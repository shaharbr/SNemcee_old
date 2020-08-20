from matplotlib import pyplot as plt
import numpy as np
import os
plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=30)  # fontsize of the figure title
plt.rcParams['font.sans-serif'] = 'Arial'



# define colormap for plotting, the colors each filter will be presented in
colormap = {'i': 'firebrick', 'r': 'tomato', 'g': 'turquoise',
            'V': 'limegreen', 'B': 'blue', 'U': 'darkorchid', 'G': 'indianred', 'R': 'tomato', 'I': 'firebrick',
            'o': 'orange', 'c': 'cyan',
            'UVW1': 'darkorchid', 'UVW2': 'darkorchid', 'UVM2': 'darkorchid'}

def add_rest_frame_days_from_discovery(SN_dict, correction_df):
    discovery_date = correction_df['discovery_date']
    z = correction_df['z']
    lightcurve = SN_dict['lightcurve']
    for source in lightcurve.keys():
        timedelta = (lightcurve[source]['df']['mjd'] - discovery_date) / (1 + z)
        lightcurve[source]['df']['t_from_discovery'] = timedelta
    SN_dict['lightcurve'] = lightcurve
    return SN_dict

# remove before-discovery data points
def remove_data_before_discovery(SN_dict):
    lightcurve = SN_dict['lightcurve']
    for source in lightcurve.keys():
        lightcurve[source]['df'] = lightcurve[source]['df'].loc[lightcurve[source]['df']['t_from_discovery'] >= 0]
    SN_dict['lightcurve'] = lightcurve
    return SN_dict



def remove_glactic_extinction(SN_dict, correction_df):
    for source in SN_dict['lightcurve'].keys():
        df = SN_dict['lightcurve'][source]['df']
        for filter in df['filter'].unique():
            df.loc[df['filter'] == filter, 'mag'] = df.loc[df['filter'] == filter, 'mag'] - correction_df[filter]
        SN_dict['lightcurve'][source]['df'] = df
    return SN_dict


def add_absolute_magnitude(SN_dict, correction_df):
    distance_modulus = correction_df['distance_modulus']
    for source in SN_dict['lightcurve'].keys():
        SN_dict['lightcurve'][source]['df']['abs_mag'] = SN_dict['lightcurve'][source]['df']['mag'] - distance_modulus
    return SN_dict

# TODO fix colormap
def lightcurve_plot(SN_dict_list, main_SN, correction_params):
    distance_modulus = correction_params.loc[main_SN]['distance_modulus']
    fig, ax = plt.subplots(1, figsize=(9, 6))
    ax2 = ax.twinx()
    names= []
    for SN in SN_dict_list:
        name = SN['Name']
        names.append(name)
        lightcurve = SN['lightcurve']
        for source in lightcurve.keys():
            df = lightcurve[source]['df']
            marker = lightcurve[source]['marker']
            linestyle = lightcurve[source]['Linestyle']
            for filter in df['filter'].unique():
                df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='abs_mag', ax=ax2,
                                                    marker=marker, linestyle=linestyle,
                                                    color=colormap[filter],
                                                    markeredgecolor='k', markeredgewidth =0.3,
                                                    label=source + ' (' + filter + ')',)
                if SN['Name'] == main_SN:
                    df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='mag', yerr='dmag', ax=ax,
                                                        marker=marker, linestyle=linestyle,
                                                        color=colormap[filter],
                                                        markeredgecolor='k', markeredgewidth=0.3)

    ax.set_title('Light-curve over time - '+str([SN['Name'] for SN in SN_dict_list]))
    ax.set_xlabel('Rest-frame days from discovery')
    ax.set_ylabel('Apparent Magnitude')
    ax.set_ylim(15, 23)
    ax2.set_ylabel('Absolute Magnitude')
    # TODO remember that the distance module difference between the y axes is hardcoded here -
    # TODO need to find way to me this automatic
    ax2.set_ylim(15 - distance_modulus, 23 - distance_modulus)
    ax2.legend(ncol=2)
    ax.invert_yaxis()
    ax.get_legend().remove()
    ax2.invert_yaxis()
    ax.tick_params(axis='both', which='major')
    ax2.tick_params(axis='both', which='major')

    fig.savefig(os.path.join('figures', 'light-curve over time ' + str(names) + '.png'))
    fig.savefig(os.path.join('figures', 'light-curve over time ' + str(names) + '.svg'))







def lightcurve_plot_shift(SN_dict, correction_params):
    distance_modulus = correction_params.loc[SN_dict['Name']]['distance_modulus']
    fig, ax = plt.subplots(1, figsize=(12, 8))
    ax2 = ax.twinx()
    lightcurve = SN_dict['lightcurve']
    shifts = {'B': '+3', 'c': '+2', 'g': '+1', 'V': '+0', 'o': '-1', 'r': '-2', 'G': '-3', 'i': '-4'}
    for filter in shifts.keys():
        for source in lightcurve.keys():
            if filter in lightcurve[source]['df']['filter'].unique():
                df = lightcurve[source]['df']
                marker = lightcurve[source]['marker']
                linestyle = lightcurve[source]['Linestyle']
                # TODO sort the legend but sorting the order in which they are plotted -only by order of filter, not depending on which source plotted first
                filter_df = df.loc[df['filter'] == filter]
                x = filter_df['t_from_discovery']
                y_abs_mag = filter_df['abs_mag'] + int(shifts[filter])

                if filter == 'V':
                    vmax = np.min(y_abs_mag)
                    # print(vmax)
                    ax2.axhline(y=vmax, color='limegreen', alpha=0.5, linestyle='--')

                ax2.errorbar(x, y_abs_mag, marker=marker, linestyle=linestyle,
                                        color=colormap[filter],
                                        markeredgecolor='k', markeredgewidth =0.3,
                                        label=filter + ' ' + str(shifts[filter]) + ' (' + source + ')',)
                y_mag = filter_df['mag'] + int(shifts[filter])
                y_mag_err = filter_df['dmag']
                ax.errorbar(x, y_mag, yerr=y_mag_err, marker=marker, linestyle=linestyle,
                         color=colormap[filter],
                         markeredgecolor='k', markeredgewidth=0.3,
                         label=filter + ' ' + str(shifts[filter]) + ' (' + source + ')', )


    ax.set_title('Light-curve over time - '+str(SN_dict['Name']) + '\n')
    ax.set_xlabel('Rest-frame days from discovery')
    ax.set_ylabel('Apparent Magnitude')
    ax.set_ylim(12, 21)
    ax2.set_ylabel('Absolute Magnitude')
    # TODO remember that the distance module difference between the y axes is hardcoded here -
    # TODO need to find way to me this automatic
    ax2.set_yticks(np.arange(int(12 - distance_modulus),  int(21 - distance_modulus + 1), 1))
    ax2.set_ylim(12 - distance_modulus, 21 - distance_modulus)
    ax2.legend(ncol=2, loc='upper right')
    ax.invert_yaxis()
    ax2.invert_yaxis()
    ax.tick_params(axis='both', which='major')
    ax2.tick_params(axis='both', which='major')

    fig.savefig(os.path.join('figures', 'light-curve over time stacked ' + str(SN_dict['Name']) + '.png'))
    fig.savefig(os.path.join('figures', 'light-curve over time stacked ' + str(SN_dict['Name']) + '.svg'))

def remove_LCO_outlier(SN_dict):
    df = SN_dict['lightcurve']['Las Cumbres']['df']
    df = df.loc[df['t_from_discovery'] < 190]
    return df


def BV_plot(SN_dict):
    df = SN_dict['lightcurve']['Las Cumbres']['df']
    B = df.loc[df['filter'] == 'B', ('t_from_discovery', 'abs_mag')]
    B['t_from_discovery'] = B['t_from_discovery'].astype(int)
    B = B.groupby('t_from_discovery').mean()
    # B = B.set_index('t_from_discovery')
    V = df.loc[df['filter'] == 'V', ('t_from_discovery', 'abs_mag')]
    V['t_from_discovery'] = V['t_from_discovery'].astype(int)
    V = V.groupby('t_from_discovery').mean()
    # V = V.set_index('t_from_discovery')
    B_V = V['abs_mag'].rsub(B['abs_mag'])
    plt.figure()
    plt.plot(B_V.index, B_V, linestyle='None', marker='o')
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.title('B-V color curve - SN2018hmx')
    plt.ylabel('B-V absolute magnitude')
    plt.xlabel('t from discovery')


