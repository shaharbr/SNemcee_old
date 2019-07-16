import json
from pandas.io.json import json_normalize
import pandas as pd
# import re


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
    lco_phot = pd.read_csv(path, sep='\t')
    # sort filter names
    lco_phot['filter'].replace({'ip': 'i', 'gp': 'g', 'rp': 'r'}, inplace=True)
    return lco_phot


def sne_catalog_phot(path):
    # import SNE catalog photometry dataset (CSV file)
    phot = pd.read_csv(path)
    phot.rename(columns={'time': 'mjd', 'magnitude': 'mag', 'e_magnitude': 'dmag', 'band': 'filter'}, inplace=True)
    return phot


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
    ztf_phot_new = pd.read_csv('ztf_dr1_lightcurve.txt', sep=r'\s+', header=50,
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


def make_SN_dict(SN_name, lightcurves_dict, z_dict, discovery_date_dict, distance_modulus_dict,
                 galactic_extinction_dict):
    SN_dict = {'lightcurve': lightcurves_dict[SN_name],
               'name': SN_name,
               'z': z_dict[SN_name],
               'discovery_date': discovery_date_dict[SN_name],
               'distance_modulus': distance_modulus_dict[SN_name],
               'galactic_extinction': galactic_extinction_dict[SN_name]}
    return SN_dict