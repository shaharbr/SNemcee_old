import numpy as np

def snec_interpolator(requested, surrounding_values, models_dict, data_days, extend_tail=False):
    if extend_tail is not False:
        model_days = np.linspace(0, 196+extend_tail, int(1+10*(196+extend_tail)))
    else:
        model_days = np.linspace(0, 196, 1961)
    params = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix']
    param_dict = {}
    for i in range(len(params)):
        param = params[i]
        below = surrounding_values[param][0]
        above = surrounding_values[param][1]
        weight_below = (above - requested[i]) / (above - below)
        param_dict[param] = {'requested': requested[i],
                                 'below': below, 'above': above,
                                 'weight_below': weight_below,
                                 'weight_above': 1 - weight_below}
    # hierarchy of nested dict: 'Mzams', 'Ni', 'E', 'R', 'K', 'Mix. take note order in filename is different
    snec_dict = {'below': {}, 'above': {}, 'requested': {}}
    for Mdir in ['below', 'above', 'requested']:
        snec_dict[Mdir] = {'below': {}, 'above': {}, 'requested': {}}
        for Nidir in ['below', 'above', 'requested']:
            snec_dict[Mdir][Nidir] = {'below': {}, 'above': {}, 'requested': {}}
            for Edir in ['below', 'above', 'requested']:
                snec_dict[Mdir][Nidir][Edir] = {'below': {}, 'above': {}, 'requested': {}}
                for Rdir in ['below', 'above', 'requested']:
                    snec_dict[Mdir][Nidir][Edir][Rdir] = {'below': {}, 'above': {}, 'requested': {}}
                    for Kdir in ['below', 'above', 'requested']:
                        snec_dict[Mdir][Nidir][Edir][Rdir][Kdir] = {'below': {}, 'above': {}, 'requested': {}}
                        for Mixdir in ['below', 'above', 'requested']:
                            snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir] = {}

    # fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    for Mdir in ['below', 'above']:
        for Nidir in ['below', 'above']:
            for Edir in ['below', 'above']:
                for Rdir in ['below', 'above']:
                    for Kdir in ['below', 'above']:
                        for Mixdir in ['below', 'above']:
                            M = param_dict['Mzams'][Mdir]
                            Ni = param_dict['Ni'][Nidir]
                            E = param_dict['E'][Edir]
                            R = param_dict['R'][Rdir]
                            K = param_dict['K'][Kdir]
                            Mix = param_dict['Mix'][Mixdir]
                            snec_model = models_dict[M][Ni][E][R][K][Mix]
                            snec_dict[Mdir][Nidir][Edir][Rdir][Kdir][Mixdir]['snec'] = \
                                np.interp(data_days, model_days, snec_model)
                        Mix_below = snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['below']['snec']
                        Mix_above = snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['above']['snec']
                        Mix_requested = Mix_below * param_dict['Mix']['weight_below'] + Mix_above * param_dict['Mix'][
                            'weight_above']
                        snec_dict[Mdir][Nidir][Edir][Rdir][Kdir]['requested']['snec'] = Mix_requested
                    K_below = snec_dict[Mdir][Nidir][Edir][Rdir]['below']['requested']['snec']
                    K_above = snec_dict[Mdir][Nidir][Edir][Rdir]['above']['requested']['snec']
                    K_requested = K_below * param_dict['K']['weight_below'] + K_above * param_dict['K'][
                        'weight_above']
                    snec_dict[Mdir][Nidir][Edir][Rdir]['requested']['requested']['snec'] = K_requested
                    # ax.plot(K_below, label='K_below')
                    # ax.plot(K_above, label='K_above')
                    # ax.plot(K_requested, label='K_requested')
                R_below = snec_dict[Mdir][Nidir][Edir]['below']['requested']['requested']['snec']
                R_above = snec_dict[Mdir][Nidir][Edir]['above']['requested']['requested']['snec']
                R_requested = R_below * param_dict['R']['weight_below'] + R_above * param_dict['R'][
                    'weight_above']
                snec_dict[Mdir][Nidir][Edir]['requested']['requested']['requested']['snec'] = R_requested
                # ax.plot(R_below, label='R_below')
                # ax.plot(R_above, label='R_above')
                # ax.plot(R_requested, label='R_requested')
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
    M_requested = M_below * param_dict['Mzams']['weight_below'] + M_above * param_dict['Mzams']['weight_above']
    snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec'] = M_requested
    # ax.plot(M_below, label='M_below')
    # ax.plot(M_above, label='M_above')
    # ax.plot(M_requested, label='M_requested')
    # ax.legend()
    # ax.set_yscale('log')
    # ax.set_ylim(2.8e+41, 1.3e+43)
    return snec_dict['requested']['requested']['requested']['requested']['requested']['requested']['snec']


