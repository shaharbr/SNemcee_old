import numpy as np

def snec_interpolator(requested, surrounding_values, models_dict, data_days, extend_tail=False):
    if extend_tail is not False:
        model_days = np.linspace(0, 200+extend_tail, int(1+10*(200+extend_tail)))
    else:
        model_days = np.linspace(0, 200, 2001)
    params = ['Mzams', 'Ni', 'E', 'Mix']
    param_dict = {}
    for i in range(len(params)):
        param = params[i]
        below = surrounding_values[param][0]
        above = surrounding_values[param][1]
        if above == below:
            weight_below = 0
        else:
            weight_below = (above - requested[i]) / (above - below)
        param_dict[param] = {'requested': requested[i],
                                 'below': below, 'above': above,
                                 'weight_below': weight_below,
                                 'weight_above': 1 - weight_below}
    # hierarchy of nested dict: 'Mzams', 'Ni', 'E', 'Mix. take note order in filename is different
    snec_dict = {'below': {}, 'above': {}, 'requested': {}}
    for Mdir in ['below', 'above', 'requested']:
        snec_dict[Mdir] = {'below': {}, 'above': {}, 'requested': {}}
        for Nidir in ['below', 'above', 'requested']:
            snec_dict[Mdir][Nidir] = {'below': {}, 'above': {}, 'requested': {}}
            for Edir in ['below', 'above', 'requested']:
                snec_dict[Mdir][Nidir][Edir] = {'below': {}, 'above': {}, 'requested': {}}
                for Mixdir in ['below', 'above', 'requested']:
                    snec_dict[Mdir][Nidir][Edir][Mixdir] = {}

    for Mdir in ['below', 'above']:
        for Nidir in ['below', 'above']:
            for Edir in ['below', 'above']:
                for Mixdir in ['below', 'above']:
                    M = param_dict['Mzams'][Mdir]
                    Ni = param_dict['Ni'][Nidir]
                    E = param_dict['E'][Edir]
                    Mix = param_dict['Mix'][Mixdir]
                    snec_model = models_dict[M][Ni][E][Mix]
                    if isinstance(snec_model, str):
                        return 'failed SN'
                    else:
                        snec_dict[Mdir][Nidir][Edir][Mixdir]['snec'] = \
                            np.interp(data_days, model_days, snec_model)
                Mix_below = snec_dict[Mdir][Nidir][Edir]['below']['snec']
                Mix_above = snec_dict[Mdir][Nidir][Edir]['above']['snec']
                Mix_requested = Mix_below * param_dict['Mix']['weight_below'] + Mix_above * param_dict['Mix'][
                    'weight_above']
                snec_dict[Mdir][Nidir][Edir]['requested']['snec'] = Mix_requested
            E_below = snec_dict[Mdir][Nidir]['below']['requested']['snec']
            E_above = snec_dict[Mdir][Nidir]['above']['requested']['snec']
            E_requested = E_below * param_dict['E']['weight_below'] + E_above * param_dict['E']['weight_above']
            snec_dict[Mdir][Nidir]['requested']['requested']['snec'] = E_requested
        Ni_below = snec_dict[Mdir]['below']['requested']['requested']['snec']
        Ni_above = snec_dict[Mdir]['above']['requested']['requested']['snec']
        Ni_requested = Ni_below * param_dict['Ni']['weight_below'] + Ni_above * param_dict['Ni']['weight_above']
        snec_dict[Mdir]['requested']['requested']['requested']['snec'] = Ni_requested
    M_below = snec_dict['below']['requested']['requested']['requested']['snec']
    M_above = snec_dict['above']['requested']['requested']['requested']['snec']
    M_requested = M_below * param_dict['Mzams']['weight_below'] + M_above * param_dict['Mzams']['weight_above']
    snec_dict['requested']['requested']['requested']['requested']['snec'] = M_requested
    return snec_dict['requested']['requested']['requested']['requested']['snec']


