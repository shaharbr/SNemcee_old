import pandas as pd
import copy
import numpy as np
import os
import re
import matplotlib.pyplot as plt

data_dir = os.path.join('..', 'all_mag_data')



def snec_interpolator(requested_list, sampled_list, data_days, data_filters):
    params = ['M', 'Ni', 'E', 'R', 'K']
    param_dict = {params[i]:{'requested': requested_list[i], 'sampled': sampled_list[i],
                             'below': 0, 'above': 0, 'weight_below': 0, 'weight_above': 0} for i in range(len(requested_list))}
    for param in param_dict.keys():
        sampled = param_dict[param]['sampled']
        requested = param_dict[param]['requested']
        print(sampled, requested)
        if requested > sampled[-1]:
            above = sampled[-1]
            print(sampled)
        else:
            above = min([sampled[i] for i in range(len(sampled)) if sampled[i] >= requested])
        if requested < sampled[0]:
            below = sampled[0]
            print(sampled)
        else:
            below = max([sampled[i] for i in range(len(sampled)) if sampled[i] <= requested])
        if (above - below) > 0:
            weight_below = (above - requested) / (above - below)
        else:
            weight_below = 1
        weight_above = 1 - weight_below
        # print(above, below, requested)
        # print(weight_below)
        # print(weight_above)
        print(below, above)
        param_dict[param]['below'] = below
        param_dict[param]['above'] = above
        param_dict[param]['weight_below'] = weight_below
        param_dict[param]['weight_above'] = weight_above


    # hierarchy of nested dict: 'M', 'Ni', 'E', 'R', 'K'
    # param_keys = list(param_dict.keys())
    param_keys = ['M', 'Ni', 'E', 'R', 'K']
    num_param = len(param_keys)
    snec_dict = {'below': {}, 'above': {}, 'requested': {}}
    for n in range(num_param - 1):
        single_snec_dict = copy.deepcopy(snec_dict)
        for i in list(single_snec_dict.keys()):
            snec_dict[i] = copy.deepcopy(single_snec_dict)

    # fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    for Mdir in ['below', 'above']:
        for Nidir in ['below', 'above']:
            for Edir in ['below', 'above']:
                for Rdir in ['below', 'above']:
                    for Kdir in ['below', 'above']:
                        list_dir = [Mdir, Nidir, Edir, Rdir, Kdir]
                        name = ''.join([param_keys[i] + str(param_dict[param_keys[i]][list_dir[i]]) + '_' for i in
                                        range(num_param)])
                        name = re.sub('_R', '_Mix3.0_R', name)
                        # remove last '_' from name
                        # TODO if something isnt working here it might be because i need to replace the
                        # TODO instantaneous objects like snec_model, K_below etc with deepcopies
                        name = name[:-1]
                        snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['name'] = name
                        mag_file = pd.read_csv(os.path.join(data_dir, name, 'magnitudes.dat'),
                                 names=['time', 'Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'], sep=r'\s+')
                        mag_file = mag_file.abs()
                        time_col = mag_file['time'] / 86400
                        snec_model_dict = {}
                        for filter in data_filters:
                            snec_model_dict[filter] = np.interp(data_days, time_col, mag_file[filter])
                        snec_model_dict['time'] = data_days
                        snec_model = pd.DataFrame(snec_model_dict)
                        snec_model = snec_model.sort_values('time')
                        snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['snec'] = snec_model
                    K_below = snec_dict[Mdir][Nidir][Edir][Rdir]['below']['snec'] # TODO do i need to specifi or do a for all filter? allso remember to add time column
                    K_above = snec_dict[Mdir][Nidir][Edir][Rdir]['above']['snec']
                    K_requested = K_below * param_dict['K']['weight_below'] + K_above * param_dict['K'][
                        'weight_above']
                    snec_dict[Mdir][Nidir][Edir][Rdir]['requested']['snec'] = K_requested
                    # ax.plot(K_below['V'], label='K_below')
                    # ax.plot(K_above['V'], label='K_above')
                    # ax.plot(K_requested['V'], label='K_requested')
                R_below = snec_dict[Mdir][Nidir][Edir]['below']['requested']['snec']
                R_above = snec_dict[Mdir][Nidir][Edir]['above']['requested']['snec']
                R_requested = R_below * param_dict['R']['weight_below'] + R_above * param_dict['R'][
                    'weight_above']
                snec_dict[Mdir][Nidir][Edir]['requested']['requested']['snec'] = R_requested
                # ax.plot(R_below, label='R_below')
                # ax.plot(R_above, label='R_above')
                # ax.plot(R_requested, label='R_requested')
            E_below = snec_dict[Mdir][Nidir]['below']['requested']['requested']['snec']
            E_above = snec_dict[Mdir][Nidir]['above']['requested']['requested']['snec']
            E_requested = E_below * param_dict['E']['weight_below'] + E_above * param_dict['E']['weight_above']
            snec_dict[Mdir][Nidir]['requested']['requested']['requested']['snec'] = E_requested
            # ax.plot(E_below, label='E_below')
            # ax.plot(E_above, label='E_above')
            # ax.plot(E_requested, label='E_requested')
        Ni_below = snec_dict[Mdir]['below']['requested']['requested']['requested']['snec']
        Ni_above = snec_dict[Mdir]['above']['requested']['requested']['requested']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['requested']['requested']['requested']['snec'] = Ni_requested
        # ax.plot(Ni_below, label='Ni_below')
        # ax.plot(Ni_above, label='Ni_above')
        # ax.plot(Ni_requested, label='Ni_requested')
    M_below = snec_dict['below']['requested']['requested']['requested']['requested']['snec']
    M_above = snec_dict['above']['requested']['requested']['requested']['requested']['snec']
    M_requested = M_below * param_dict['M']['weight_below'] + M_above * param_dict['M']['weight_above']
    snec_dict['requested']['requested']['requested']['requested']['requested']['snec'] = M_requested
    # ax.plot(M_below, label='M_below')
    # ax.plot(M_above, label='M_above')
    # ax.plot(M_requested, label='M_requested')
    # ax.legend()
    # ax.set_yscale('log')
    # ax.set_ylim(2.8e+41, 1.3e+43)
    return snec_dict['requested']['requested']['requested']['requested']['requested']['snec']


# requested = [13.5, 0.15, 1.8, 3, 1000, 50]
# sampled = [Mzams_range, Ni_range, E_final_range, Mix_range, R_range, K_range]

# interp_model = snec_interpolator(requested, sampled)


# print(interp_model)

plt.show()

