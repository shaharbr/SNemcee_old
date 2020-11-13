import json
from pandas.io.json import json_normalize
import pandas as pd
import re
import os
import numpy as np

# from scipy.optimize import curve_fit


def convert_to_mjd(datetime_series, from_datetime=False):
    # if origin format is datetime, convert first to JD
    if from_datetime:
        datetime_series = pd.to_datetime(datetime_series)
        if isinstance(datetime_series, pd.Series):
            for i in range(len(datetime_series)):
                datetime_series[i] = datetime_series[i].to_julian_date()
        else:
            datetime_series = datetime_series.to_julian_date()
    # convert from JD to MJD
    datetime_series = datetime_series - 2400000.5
    return datetime_series


def lco_phot(path):
    # import LCO photometry dataset (TSV file)
    lco_phot = pd.read_csv(path, sep=r'\s+')
    lco_phot['mjd'] = convert_to_mjd(lco_phot['mjd'])
    # sort filter names
    # lco_phot['filter'].replace({'ip': 'i', 'gp': 'g', 'rp': 'r'}, inplace=True)
    return lco_phot


def sn18aad_phot(path):
    # import LCO photometry dataset (TSV file)
    phot = pd.read_csv(path, sep=r'\s+', names=['date', 'mjd', 'mag', 'dmag', 'telescope', 'filter'])
    #TODO sort the title here
    phot['mjd'] = convert_to_mjd(phot['mjd'])
    return phot


def ksp_phot(path):
    # import KSP photometry dataset (TSV file)
    ksp_phot = pd.read_csv(path, sep=r'\s+', names=['mjd', 'filter', 'mag', 'dmag'])
    return ksp_phot


def sne_catalog_phot(path):
    # import SNE catalog photometry dataset (CSV file)
    phot = pd.read_csv(path)
    phot.rename(columns={'time': 'mjd', 'magnitude': 'mag', 'e_magnitude': 'dmag', 'band': 'filter'}, inplace=True)
    return phot


def atlas_flux(path, filter):
    df = pd.read_csv(path, sep='\t', header=0)
    df['filter'] = filter
    df = df.loc[df['flux'] > (3 * df['flux_error'])]
    df['mag'] = -2.5 * np.log10(df['flux'] * 10**(-6)) + 8.9
    df['dmag'] = -2.5 * np.log10(np.exp(1)) / df['flux'] * df['flux_error']
    # filter nondetection
    return df


def gaia_phot(path):
    # import Gaia photometry dataset (CSV file)
    gaia_phot = pd.read_csv(path, header=1)
    # standardize column names
    gaia_phot.rename(columns={'JD': 'mjd', 'averagemag.': 'mag'}, inplace=True)
    # convert JD to MJD
    gaia_phot['mjd'] = convert_to_mjd(gaia_phot['mjd'])
    # fill missing columns
    gaia_phot['dmag'] = 0
    # sort filter names
    gaia_phot['filter'] = 'G'
    return gaia_phot


def atlas_phot(path):
    # import ATLAS photometry dataset (TSV file)
    atlas_phot = pd.read_csv(path, sep='\t')
    # standardize column names
    atlas_phot.rename(columns={'Obs-date': 'mjd', 'Mag. / Flux': 'mag', 'Err': 'dmag'}, inplace=True)
    atlas_phot['mjd'] = convert_to_mjd(atlas_phot['mjd'], from_datetime=True)
    # sort filter names
    atlas_phot['filter'] = 'o'
    return atlas_phot


def ztf_phot(path):
    # import ZTF photometry dataset (JSON file)
    file = open(path)
    file = json.load(file)
    ztf_phot = json_normalize(file['results'])[
        ['candidate.jd', 'candidate.magpsf', 'candidate.sigmapsf', 'candidate.filter']]
    # standardize column names
    ztf_phot.rename(columns={'candidate.jd': 'mjd', 'candidate.magpsf': 'mag', 'candidate.sigmapsf': 'dmag',
                             'candidate.filter': 'filter'}, inplace=True)
    # convert JD to MJD
    ztf_phot['mjd'] = convert_to_mjd(ztf_phot['mjd'])
    return ztf_phot


def ztf_phot_new(path):
    # import re-analyzed ZTF photometry dataset (ASCII file)
    ztf_phot_new = pd.read_csv(os.path.join('data', 'ztf_dr1_lightcurve.txt'), sep=r'\s+', header=50,
                               usecols=['mjd|', 'mag|', 'hjd|', 'catflags|'])
    ztf_phot_new.drop([0, 1, 2], inplace=True)
    # standardize column names and fill missing fields
    ztf_phot_new.rename(columns={'hjd|': 'mjd', 'mjd|': 'mag', 'mag|': 'dmag', 'catflags|': 'filter'},
                        inplace=True)
    # sort out types
    ztf_phot_new['mjd'], ztf_phot_new['mag'], ztf_phot_new['dmag'], ztf_phot_new['filter'] = \
        ztf_phot_new['mjd'].astype('float64'), ztf_phot_new['mag'].astype('float64'), \
        ztf_phot_new['dmag'].astype('float64'), ztf_phot_new['filter'].astype('str'),
    # sort filter names
    ztf_phot_new['filter'].replace({'zg': 'g', 'zr': 'r'}, inplace=True)
    return ztf_phot_new


def leonard_phot(path):
    # import photometry dataset from Leonard et al 2002 (TSV file)
    leonard_phot = pd.read_csv(path, sep='\s+', header=None,
                               names=['t_from_discovery', 'U', 'B', 'V', 'R', 'I'], parse_dates=['t_from_discovery'])
    # restructure leonard data as the others - all mag in one column, and another column which indicates
    #  which filter its from
    leonard_phot = leonard_phot.melt(id_vars=['t_from_discovery'], var_name='filter', value_name='mag')
    leonard_phot['t_from_discovery'] = leonard_phot['t_from_discovery'].astype('float64')
    leonard_phot = leonard_phot.loc[leonard_phot['t_from_discovery'] < 300]
    leonard_phot = leonard_phot.loc[leonard_phot['mag'] < 90]
    leonard_phot['dmag'] = 0
    return leonard_phot

def date_to_time_from_discovery(date, discovery_date):
    timedelta = (date - discovery_date)
    return timedelta

def t_to_restframe(t, z):
    t = t / (1 + z)
    return t





def remove_last_days(df, day_limit):
    df = df.loc[df['t_from_discovery'] < day_limit]
    return df



def SN14hls_expans_v(path):
    # import expansion velocity file for iPTF14hls
    expans_v_df = pd.read_csv(path, header=0)
    expans_v_df['JD'] = pd.to_datetime(expans_v_df['JD'], unit='D', origin='julian')
    expans_v_df.rename(columns={'JD': 'datetime', 'Velocity [km/s]': 'absorption_mean_velocity', 'Line': 'line',
                                       'Velocity_Error [km/s]': 'absorption_std_velocity'}, inplace=True)
    return expans_v_df

# def fit_with_errors(fit_function, x, y, err=False):
#     popt, pcov = curve_fit(fit_function, x, y, sigma=err, absolute_sigma=True)
#     return popt, pcov


def dict_to_csv(dict, filename):
    pd.DataFrame(dict, index=[0]).to_csv(os.path.join('data', filename))

def correction_dict_to_csv(z_dict, discovery_date_dict, distance_modulus_dict, galactic_extinction_dict):
    for SN in z_dict.keys():
        dict = {'Name': SN}
        dict['z'] = z_dict[SN]
        dict['discovery_date'] = convert_to_mjd(discovery_date_dict[SN], from_datetime=True)
        dict['dm'] = distance_modulus_dict[SN]
        for filter in galactic_extinction_dict[SN].keys():
            dict[filter] = galactic_extinction_dict[SN][filter]
        pd.DataFrame(dict, index=[0]).to_csv(os.path.join('data', SN + '_correction.csv'))




def make_alllightcurve_df(SN_dict):
    sources = SN_dict['lightcurve'].keys()
    all_df = []
    for source in sources:
        df = SN_dict['lightcurve'][source]['df'].copy(deep=True)
        df = df.loc[df['mag'] > 0]
        df['source'] = source
        all_df.append(df)
    all_df = pd.concat(all_df, sort=False)
    all_df = all_df[['mjd', 'mag', 'dmag', 'filter', 'source']]
    all_df.reset_index(inplace=True, drop=True)
    return all_df


def save_ascii(dataframe, filename):
    df = dataframe.copy(deep=True)
    df.rename(columns={'mjd': 'MJD', 'filter': 'filt'}, inplace=True)
    df = df.loc[((df['source'] == 'Las Cumbres') |
                 (df['source'] == 'P60') |
                 (df['source'] == 'Keck') |
                 (df['source'] == 'Arizona')|
                 (df['source'] == 'SNe catalog')) &
                ((df['filt'] == 'B') |
                 (df['filt'] == 'V') |
                 (df['filt'] == 'R') |
                 (df['filt'] == 'g') |
                 (df['filt'] == 'r') |
                 (df['filt'] == 'i'))]
    df['nondet |'] = 'FALSE |'
    header_row = pd.DataFrame(np.array(df.keys()).reshape(1, 6), columns=list(df.keys()))
    df = pd.concat([header_row, df], ignore_index=True)
    def equalspace(x):
        x = str(x)
        len_x = len(x)
        x = '|' + ' '*(30-len_x) + x
        return x
    for column in df.keys():
        df[column] = df[column].apply(equalspace)
    df['merged'] = df[df.columns[0:]].apply(lambda x: ''.join(x),axis=1)
    df['merged'].to_csv(filename, index=False, header=False)
    return df
