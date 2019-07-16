from matplotlib import pyplot as plt


def add_rest_frame_days_from_discovery(SN_dict):
    discovery_date = SN_dict['discovery_date']
    z = SN_dict['z']
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


def remove_glactic_extinction(SN_dict):
    lightcurve = SN_dict['lightcurve']
    galactic_extinction_values = SN_dict['galactic_extinction']
    for source in lightcurve.keys():
        df = lightcurve[source]['df']
        for filter in df['filter'].unique():
            df.loc[df['filter'] == filter, 'mag'] = df.loc[df['filter'] == filter, 'mag'] - galactic_extinction_values[filter]
    SN_dict['lightcurve'] = lightcurve
    return SN_dict


def add_absolute_magnitude(SN_dict):
    lightcurve = SN_dict['lightcurve']
    distance_modulus = SN_dict['distance_modulus']
    for source in lightcurve.keys():
        lightcurve[source]['df']['abs_mag'] = lightcurve[source]['df']['mag'] - distance_modulus
    SN_dict['lightcurve'] = lightcurve
    return SN_dict

# colormap
def lightcurve_plot(SN_dict_list, main_SN, V50_line=False):
    fig, ax = plt.subplots(1, figsize=(9, 6))
    ax2 = ax.twinx()
    for SN in SN_dict_list:
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
                if SN['name'] == main_SN:
                    distance_modulus = SN['distance_modulus']
                    df.loc[df['filter'] == filter].plot(x='t_from_discovery', y='mag', yerr='dmag', ax=ax,
                                                        marker=marker, linestyle=linestyle,
                                                        color=colormap[filter],
                                                        markeredgecolor='k', markeredgewidth=0.3)

    ax.set_title('Light-curve over time - '+str([SN['name'] for SN in SN_dict_list]), fontsize=16)
    ax.set_xlabel('Time since discovery (rest-frame days)', size=16)
    ax.set_ylabel('Apparent Magnitude', size=16)
    ax.set_ylim(14, 27)
    ax2.set_ylabel('Absolute Magnitude', size=16)
    # TODO remember that the distance module difference between the y axes is hardcoded here -
    # TODO need to find way to me this automatic
    ax2.set_ylim(14 - distance_modulus, 27 - distance_modulus)
    ax2.legend(ncol=2)
    ax.invert_yaxis()
    ax.get_legend().remove()
    ax2.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)

    # fig.savefig('light-curve over time - 2018hmx vs 1999em' + '.png')

def remove_LCO_outlier(SN_dict):
    df = SN_dict['lightcurve']['LCO']['df']
    df = df.loc[df['t_from_discovery'] < 190]
    return df

