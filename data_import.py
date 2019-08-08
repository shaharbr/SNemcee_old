import json
from pandas.io.json import json_normalize
import pandas as pd
import re
import os

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

def date_to_time_from_discovery(date, discovery_date):
    timedelta = (date - discovery_date)
    return timedelta

def t_to_restframe(t, z):
    t = t / (1 + z)
    return t


def make_SN_dict(SN_name, lightcurves_dict=False, z_dict=False, discovery_date_dict=False,
                 distance_modulus_dict=False, galactic_extinction_dict=False,
                 spectra=False, expansion_velocities=False):

    # organize in SN_dict
    fields = [lightcurves_dict, z_dict, discovery_date_dict,
                 distance_modulus_dict, galactic_extinction_dict,
                 spectra, expansion_velocities]
    keys = ['lightcurve', 'z','discovery_date', 'distance_modulus', 'galactic_extinction',
               'spectra', 'expansion_velocities']

    # convert date to MJD
    discovery_date_dict[SN_name] = convert_to_mjd(discovery_date_dict[SN_name], from_datetime=True)
    SN_dict = {}
    SN_dict['name'] = SN_name
    for i in range(len(fields)):
        if fields[i]:
            SN_dict[keys[i]] = fields[i][SN_name]
        else:
            SN_dict[keys[i]] = ''
    return SN_dict


def remove_last_days(df, day_limit):
    df = df.loc[df['t_from_discovery'] < day_limit]
    return df

def LCO_spect(dir_path):
    # import all LCO spectra ascii files for 2018hmx and organize in a dictionary
    folder_join = os.path.join
    filenames = os.listdir(dir_path)
    # reading and merging
    LCO_spect_dict = {}
    for file in filenames:
        # extract date of spectrum measurment from filename
        # TODO make this re condition more general for all SNe
        date = re.sub('.*hmx|.*aad|_|-|P60.*|v1|[a-z]|[A-Z]|\..*', '', os.path.basename(file))
        # transform to standard datetime format
        date = pd.to_datetime(date)
        # convert date to MJD
        date = convert_to_mjd(date, from_datetime=True)
        # add as element in dict
        if 'ZTF' in os.path.basename(file):
            LCO_spect_dict[date] = {
                'df': pd.read_csv(folder_join(dir_path, file), sep=' ', names=["x", "y", 'dy'], header=180)}
        else:
            LCO_spect_dict[date] = {'df': pd.read_csv(folder_join(dir_path, file), sep=' ', names=["x", "y"])}
        if 'redblu' in os.path.basename(file):
            LCO_spect_dict[date]['telescope'] = 'LCO'
        elif 'ZTF' in os.path.basename(file):
            LCO_spect_dict[date]['telescope'] = 'ZTF'
        elif 'HET' in os.path.basename(file):
            LCO_spect_dict[date]['telescope'] = 'HET'
        else:
            LCO_spect_dict[date]['telescope'] = 'ND'
    return LCO_spect_dict

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