import pandas as pd
import copy
import numpy as np
import os
import matplotlib.pyplot as plt

# def getshape(d):
#     if isinstance(d, dict):
#         print(len(d.keys()))
#         first = list(d.keys())[0]
#         getshape(d[first])
#
# print(getshape(snec_dict))
# data_dir = '/home/sbracha/SNEC/all_data/'

data_dir = 'all_data'

def simple_interpolator(M_requested, Ni_requested, M_sampled, Ni_sampled):
    requested = [M_requested, Ni_requested]
    sampled = [M_sampled, Ni_sampled]
    params = ['M', 'Ni']
    param_dict = {params[i]:{'requested': requested[i], 'sampled': sampled[i],
                             'below': 0, 'above': 0, 'weight_below': 0, 'weight_above': 0} for i in range(len(requested))}
    for param in param_dict.keys():
        sampled = param_dict[param]['sampled']
        requested = param_dict[param]['requested']
        below = max([sampled[i] for i in range(len(sampled)) if sampled[i] <= requested])
        above = min([sampled[i] for i in range(len(sampled)) if sampled[i] >= requested])
        weight_below = (above - requested) / (above - below)
        weight_above = 1 - weight_below
        print(above, below, requested)
        print(weight_below)
        print(weight_above)

        param_dict[param]['below'] = below
        param_dict[param]['above'] = above
        param_dict[param]['weight_below'] = weight_below
        param_dict[param]['weight_above'] = weight_above


    # hierarchy of nested dict: 'M', 'Ni'
    param_keys = list(param_dict.keys())
    num_param = len(param_keys)
    snec_dict = {'below': {}, 'above': {}, 'requested': {}}
    for n in range(num_param - 1):
        single_snec_dict = copy.deepcopy(snec_dict)
        for i in list(single_snec_dict.keys()):
            snec_dict[i] = copy.deepcopy(single_snec_dict)

    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    for Mdir in ['below', 'above']:
        for Nidir in ['below', 'above']:
            list_dir = [Mdir, Nidir]
            name = ''.join([param_keys[i] + str(param_dict[param_keys[i]][list_dir[i]]) + '_' for i in
                            range(num_param)])
            # remove last '_' from name
            # TODO if something isnt working here it might be because i need to replace the
            # TODO instantaneous objects like snec_model, K_below etc with deepcopies
            name = name[:-1]
            snec_dict[Mdir][Nidir]['name'] = name
            snec_model = pd.read_csv(os.path.join(data_dir, name, 'lum_observed.dat'),
                                     names=['t_from_discovery', 'Lum'], sep=r'\s+')
            snec_model = snec_model['Lum']
            snec_dict[Mdir][Nidir]['snec'] = snec_model
        Ni_below = snec_dict[Mdir]['below']['snec']
        Ni_above = snec_dict[Mdir]['above']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['snec'] = Ni_requested
        ax.plot(Ni_below, label='M'+str(param_dict['M'][Mdir])+'_'+'Ni'+str(param_dict['Ni']['below']))
        ax.plot(Ni_above, label='M'+str(param_dict['M'][Mdir])+'_'+'Ni'+str(param_dict['Ni']['above']))
    ax.plot(Ni_requested, label='M'+str(param_dict['M'][Mdir])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    M_below = snec_dict['below']['requested']['snec']
    M_above = snec_dict['above']['requested']['snec']
    M_requested = M_below * param_dict['M']['weight_below'] + M_above * param_dict['M']['weight_above']
    snec_dict['requested']['requested'] = M_requested
    ax.plot(M_below, label='M'+str(param_dict['M']['below'])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    ax.plot(M_above, label='M'+str(param_dict['M']['above'])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    ax.plot(M_requested, label='M'+str(param_dict['M']['requested'])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    ax.legend()
    ax.set_yscale('log')
    ax.set_ylim(2.8e+41, 1.3e+43)


simple_interpolator(13.5, 0.13, [13.0, 15.0], [0.1, 0.2])
plt.show()

def snec_interpolator(M_requested, Ni_requested, E_requested, Mix_requested, R_requested, K_requested, M_sampled, Ni_sampled, E_sampled, Mix_sampled, R_sampled, K_sampled):
    requested = [M_requested, Ni_requested, E_requested, Mix_requested, R_requested, K_requested]
    sampled = [M_sampled, Ni_sampled, E_sampled, Mix_sampled, R_sampled, K_sampled]
    params = ['M', 'Ni', 'E', 'Mix', 'R', 'K']
    param_dict = {params[i]:{'requested': requested[i], 'sampled': sampled[i],
                             'below': 0, 'above': 0, 'weight_below': 0, 'weight_above': 0} for i in range(len(requested))}
    for param in param_dict.keys():
        sampled = param_dict[param][sampled]
        requested = param_dict[param][requested]
        below = max(sampled < requested)
        above = min(sampled > requested)
        weight_below = (above - below) / (requested - below)
        weight_above = 1 - weight_below

        param_dict[param]['below'] = below
        param_dict[param]['above'] = above
        param_dict[param]['weight_below'] = weight_below
        param_dict[param]['weight_above'] = weight_above

    # init_below = ''.join([param+str(param['below'])+'_' for param in param_dict])
    # init_above = ''.join([param + str(param['above']) + '_' for param in param_dict])
    # below_snec = pd.read_csv(data_dir + init_below + r'\lum_observed.dat', names=['t_from_discovery', 'Lum'],
    #                          sep=r'\s+')
    # above_snec = pd.read_csv(data_dir + init_above + r'\lum_observed.dat', names=['t_from_discovery', 'Lum'],
    #                          sep=r'\s+')


    # hierarchy of nested dict: 'M', 'Ni', 'E', 'Mix', 'R', 'K'
    param_keys = list(param_dict.keys())
    num_param = len(param_keys)
    snec_dict = {'below': 0, 'above': 0, 'requested': 0}
    for n in range(num_param - 1):
        single_snec_dict = copy.deepcopy(snec_dict)
        for i in list(single_snec_dict.keys()):
            snec_dict[i] = copy.deepcopy(single_snec_dict)

    print(snec_dict)

    # M
    for Mdir in ['below', 'above']:
        for Nidir in ['below', 'above']:
            for Edir in ['below', 'above']:
                for Mixdir in ['below', 'above']:
                    for Rdir in ['below', 'above']:
                        for Kdir in ['below', 'above']:
                            list_dir = [Mdir, Nidir, Edir, Mixdir, Rdir, Kdir]
                            name = ''.join([param_keys[i] + param_dict[param_keys[i][list_dir[i]]] + '_' for i in
                                            range(num_param)])
                            # remove last '_' from name
                            # TODO if something isnt working here it might be because i need to replace the
                            # TODO instantaneous objects like snec_model, K_below etc with deepcopies
                            name = name[:-1]
                            snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir][Kdir]['name'] = name
                            snec_model = pd.read_csv(data_dir + name + r'\lum_observed.dat',
                                                     names=['t_from_discovery', 'Lum'], sep=r'\s+')
                            snec_model = snec_model['Lum']
                            snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir][Kdir]['snec'] = snec_model
                        K_below = snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir]['below']['snec']
                        K_above = snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir]['above']['snec']
                        K_requested = K_below * param_dict['K']['weight_below'] + K_above * param_dict['K']['weight_above']
                        snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir]['requested']['snec'] = K_requested
                    R_below = snec_dict[Mdir][Nidir][Edir][Mixdir]['below']['requested']['snec']
                    R_above = snec_dict[Mdir][Nidir][Edir][Mixdir]['above']['requested']['snec']
                    R_requested = R_below * param_dict['R']['weight_below'] + R_above * param_dict['R']['weight_above']
                    snec_dict[Mdir][Nidir][Edir][Mixdir]['requested']['requested']['snec'] = R_requested
                Mix_below = snec_dict[Mdir][Nidir][Edir]['below']['requested']['requested']['snec']
                Mix_above = snec_dict[Mdir][Nidir][Edir]['above']['requested']['requested']['snec']
                Mix_requested = Mix_below * param_dict['Mix']['weight_below'] + Mix_above * param_dict['Mix']['weight_above']
                snec_dict[Mdir][Nidir][Edir]['requested']['requested']['requested']['snec'] = Mix_requested
            E_below = snec_dict[Mdir][Nidir]['below']['requested']['requested']['requested']['snec']
            E_above = snec_dict[Mdir][Nidir]['above']['requested']['requested']['requested']['snec']
            E_requested = E_below * param_dict['E']['weight_below'] + E_above * param_dict['E']['weight_above']
            snec_dict[Mdir][Nidir]['requested']['requested']['requested']['requested']['snec'] = E_requested
        Ni_below = snec_dict[Mdir]['below']['requested']['requested']['requested']['requested']['snec']
        Ni_above = snec_dict[Mdir]['above']['requested']['requested']['requested']['requested']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['requested']['requested']['requested']['requested']['snec'] = Ni_requested
    M_below = snec_dict['below']['requested']['requested']['requested']['requested']['requested']['snec']
    M_above = snec_dict['above']['requested']['requested']['requested']['requested']['requested']['snec']
    M_requested = M_below * param_dict['M']['weight_below'] + M_above * param_dict['M']['weight_above']
    snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec'] = M_requested
