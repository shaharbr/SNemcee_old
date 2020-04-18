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

data_dir = os.path.join(r'C:\Users','Shahar','Documents', 'all_lum_data')



# Mzams_range = [13.0, 15.0]
# Ni_range = [0.1, 0.2]
# E_final_range = [1.8, 1.8]
# Mix_range = [3.0, 3.0]
# R_range = [1000, 1000]
# K_range = [50, 50]



def simple_interpolator(requested, sampled):
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
        # print(weight_below)
        # print(weight_above)

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

    # fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
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
            time_col = snec_model['t_from_discovery'] / 86400
            snec_model = snec_model['Lum']
            interp_days = np.arange(15, 170, 1)
            snec_model = np.interp(interp_days, time_col, snec_model)
            snec_dict[Mdir][Nidir]['snec'] = snec_model
        Ni_below = snec_dict[Mdir]['below']['snec']
        Ni_above = snec_dict[Mdir]['above']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['snec'] = Ni_requested
        # ax.plot(Ni_below, label='M'+str(param_dict['M'][Mdir])+'_'+'Ni'+str(param_dict['Ni']['below']))
        # ax.plot(Ni_above, label='M'+str(param_dict['M'][Mdir])+'_'+'Ni'+str(param_dict['Ni']['above']))
    # ax.plot(Ni_requested, label='M'+str(param_dict['M'][Mdir])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    M_below = snec_dict['below']['requested']['snec']
    M_above = snec_dict['above']['requested']['snec']
    M_requested = M_below * param_dict['M']['weight_below'] + M_above * param_dict['M']['weight_above']
    snec_dict['requested']['requested'] = M_requested
    # print(type(snec_dict['requested']['requested']))
    # time_col = snec_model['t_from_discovery']
    # ax.plot(M_below, label='M'+str(param_dict['M']['below'])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    # ax.plot(M_above, label='M'+str(param_dict['M']['above'])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    # ax.plot(M_requested, label='M'+str(param_dict['M']['requested'])+'_'+'Ni'+str(param_dict['Ni']['requested']))
    # ax.legend()
    # ax.set_yscale('log')
    # ax.set_ylim(2.8e+41, 1.3e+43)
    return snec_dict['requested']['requested']






def snec_interpolator(requested_list, sampled_list):
    params = ['M', 'Ni', 'E', 'Mix', 'R', 'K']
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


    # hierarchy of nested dict: 'M', 'Ni', 'E', 'Mix', 'R', 'K'
    param_keys = list(param_dict.keys())
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
                for Mixdir in ['below', 'above']:
                    for Rdir in ['below', 'above']:
                        for Kdir in ['below', 'above']:
                            list_dir = [Mdir, Nidir, Edir, Mixdir, Rdir, Kdir]
                            name = ''.join([param_keys[i] + str(param_dict[param_keys[i]][list_dir[i]]) + '_' for i in
                                            range(num_param)])


                            # remove last '_' from name
                            # TODO if something isnt working here it might be because i need to replace the
                            # TODO instantaneous objects like snec_model, K_below etc with deepcopies
                            name = name[:-1]

                            snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir][Kdir]['name'] = name

                            # snec_model = pd.read_csv(r'all_data\M13.0_Ni0.1_E1.8_Mix3.0_R1000_K50\lum_observed.dat')

                            snec_model = pd.read_csv(os.path.join(data_dir, name, 'lum_observed.dat'),
                                     names=['t_from_discovery', 'Lum'], sep=r'\s+')

                            time_col = snec_model['t_from_discovery'] / 86400
                            snec_model = snec_model['Lum']
                            interp_days = np.arange(15, 170, 1)
                            snec_model = np.interp(interp_days, time_col, snec_model)
                            snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir][Kdir]['snec'] = snec_model
                        K_below = snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir]['below']['snec']
                        K_above = snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir]['above']['snec']
                        K_requested = K_below * param_dict['K']['weight_below'] + K_above * param_dict['K'][
                            'weight_above']
                        snec_dict[Mdir][Nidir][Edir][Mixdir][Rdir]['requested']['snec'] = K_requested
                        # ax.plot(K_below, label='K_below')
                        # ax.plot(K_above, label='K_above')
                        # ax.plot(K_requested, label='K_requested')
                    R_below = snec_dict[Mdir][Nidir][Edir][Mixdir]['below']['requested']['snec']
                    R_above = snec_dict[Mdir][Nidir][Edir][Mixdir]['above']['requested']['snec']
                    R_requested = R_below * param_dict['R']['weight_below'] + R_above * param_dict['R'][
                        'weight_above']
                    snec_dict[Mdir][Nidir][Edir][Mixdir]['requested']['requested']['snec'] = R_requested
                    # ax.plot(R_below, label='R_below')
                    # ax.plot(R_above, label='R_above')
                    # ax.plot(R_requested, label='R_requested')
                Mix_below = snec_dict[Mdir][Nidir][Edir]['below']['requested']['requested']['snec']
                Mix_above = snec_dict[Mdir][Nidir][Edir]['above']['requested']['requested']['snec']
                Mix_requested = Mix_below * param_dict['Mix']['weight_below'] + Mix_above * param_dict['Mix'][
                    'weight_above']
                snec_dict[Mdir][Nidir][Edir]['requested']['requested']['requested']['snec'] = Mix_requested
                # ax.plot(Mix_below, label='Mix_below')
                # ax.plot(Mix_above, label='Mix_above')
                # ax.plot(Mix_requested, label='Mix_requested')
            E_below = snec_dict[Mdir][Nidir]['below']['requested']['requested']['requested']['snec']
            E_above = snec_dict[Mdir][Nidir]['above']['requested']['requested']['requested']['snec']
            E_requested = E_below * param_dict['E']['weight_below'] + E_above * param_dict['E']['weight_above']
            snec_dict[Mdir][Nidir]['requested']['requested']['requested']['requested']['snec'] = E_requested
            # ax.plot(E_below, label='E_below')
            # ax.plot(E_above, label='E_above')
            # ax.plot(E_requested, label='E_requested')
        Ni_below = snec_dict[Mdir]['below']['requested']['requested']['requested']['requested']['snec']
        Ni_above = snec_dict[Mdir]['above']['requested']['requested']['requested']['requested']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['requested']['requested']['requested']['requested']['snec'] = Ni_requested
        # ax.plot(Ni_below, label='Ni_below')
        # ax.plot(Ni_above, label='Ni_above')
        # ax.plot(Ni_requested, label='Ni_requested')
    M_below = snec_dict['below']['requested']['requested']['requested']['requested']['requested']['snec']
    M_above = snec_dict['above']['requested']['requested']['requested']['requested']['requested']['snec']
    M_requested = M_below * param_dict['M']['weight_below'] + M_above * param_dict['M']['weight_above']
    snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec'] = M_requested
    # ax.plot(M_below, label='M_below')
    # ax.plot(M_above, label='M_above')
    # ax.plot(M_requested, label='M_requested')
    # ax.legend()
    # ax.set_yscale('log')
    # ax.set_ylim(2.8e+41, 1.3e+43)
    return snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec']


# requested = [13.5, 0.15, 1.8, 3, 1000, 50]
# sampled = [Mzams_range, Ni_range, E_final_range, Mix_range, R_range, K_range]

# interp_model = snec_interpolator(requested, sampled)


# print(interp_model)

# plt.show()



def old_snec_interpolator(M_requested, Ni_requested, E_requested, Mix_requested, R_requested, K_requested, M_sampled, Ni_sampled, E_sampled, Mix_sampled, R_sampled, K_sampled):
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

                            snec_model = pd.read_csv(os.path.join(data_dir, name, 'lum_observed.dat'),
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
