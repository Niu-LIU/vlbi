#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_ts.py
"""
Created on Wed Feb  5 13:44:22 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Read the radio source coordinate time-series and return an astropy.table.Table object.
"""

from astropy.table import Table, Column
from astropy.time import Time

import astropy.units as u
import numpy as np
# My mudoles
from my_progs.catalog.read_icrf import read_icrf3
from my_progs.catalog.pos_diff import pa_calc


# -----------------------------  FUNCTIONS -----------------------------
def mjd2mjy(mjd):
    """Convert MJD to MJY.

    Parameter
    ---------
    mjd : 1darray-like
        mean Julian day

    Return
    ------
    mjy : 1darray-like
        mean Julian year
    """

    epoch = Time(mjd, format="mjd")
    mjy = epoch.jyear

    return mjy


def read_ts(datafile, datadir="../data"):
    """Read radio source coordinate time-series data

    Parameters
    ----------
    datafile: string
        data file name
    datadir: string
        directory to store the data file

    Returns
    -------
    coordts : astropy.table.Table object
        coordinate time series
    """

    coordts = Table.read("{}/{}".format(datadir, datafile), format="ascii",
                         names=["mjd", "ra", "dec", "ra_err", "dec_err", "ra_dec_corr",
                                "used_obs", "iers_name", "ivs_name", "db_name"])

    # Add unit
    coordts["ra"].unit = u.deg
    coordts["dec"].unit = u.deg
    coordts["ra_err"].unit = u.mas
    coordts["dec_err"].unit = u.mas

    mjy = mjd2mjy(coordts["mjd"])
    mjy = Column(mjy, name="mjy")
    coordts.add_column(mjy)

    return coordts


def get_sou_ts(souname, datadir="../data"):
    """Read radio source coordinate time-series data

    Parameters
    ----------
    souname: string
        source name of IERS designation or file name with the data
    datadir: string
        directory to store the data file

    Returns
    -------
    coordts : astropy.table.Table object
        coordinate time series
    """

    datafile = "{}.txt".format(souname)
    coordts = read_ts(datafile, datadir)

    return coordts


def ts_stat(coordts):
    """Simple statistics of coordinate time-series

    Parameter
    ---------
    souname: string
        source name of IERS designation or file name with the data
    """

    stat = {}

    # Number of used sessions
    numses = len(coordts)
    stat["used_ses"] = numses

    # Number of used observations
    numobs = np.sum(coordts["used_obs"])
    stat["used_obs"] = numobs

    #### TO-DO ####
    # MAYBE PRINT THESE STATISTICS
    #### TO-DO ####


def sou_ts_stat(souname):
    """Simple statistics of coordinate time-series

    Parameter
    ---------
    souname: string
        source name
    """

    coordts = get_sou_ts(souname)

    ts_stat(coordts)

    # NOT CLEAR YET ABOUT WHAT TO DO NEXT


def get_icrf3_pos(souname, wv="sx"):
    """Get radio source position from the ICRF3 catalog

    Parameters
    ----------
    souname: string
        source name
    wv: str,
        wavelength, SX, K, or XKa

    Returns
    -------
    pos : a dict object
        key: "ra", "dec"
        value: ra, dec in degree
    """

    icrf3 = read_icrf3(wv=wv)

    # Find the source position
    mask = (icrf3["iers_name"] == souname)
    sou = icrf3[mask]

    if len(sou):
        ra = sou["ra"][0]
        dec = sou["dec"][0]

    else:
        print("There is no entry for {} in ICRF3 catalog".format(souname))
        print("Using the mean position")
        return False

    pos = {"ra": ra, "dec": dec}

    return pos


def calc_ts_oft(souname, coordts):
    """Calculate the time series offset wrt. ICRF3 position

    Parameters
    ----------
    souname: string
        source name
    coordts : astropy.table.Table object
        coordinate time series

    Return
    ------
    dpos : a dict object
        key: "dra", "ddec", "rho", "pa"
        value: dra, ddec, rho in mas
               pa in deg
    """

    pos = get_icrf3_pos(souname)
    #
    if pos is False:
        ra0 = pos["ra"]
        dec0 = pos["dec"]
    else:
        ra0 = np.mean(coordts["ra"])
        dec0 = np.mean(coordts["dec"])

    dra = coordts["ra"] - ra0
    ddec = coordts["dec"] - dec0

    dra.unit = u.deg
    ddec.unit = u.deg
    dra.convert_unit_to(u.mas)
    ddec.convert_unit_to(u.mas)

    rho = np.sqrt(dra**2 + ddec**2)
#     pax, pay = pa_calc(np.array(dra), np.array(ddec))
    pa = pa_calc(np.array(dra), np.array(ddec))

    dpos = {"dra": dra, "ddec": ddec, "rho": rho, "pa": pa}

    return dpos


def get_ts(souname, datadir="../data"):
    """get the time series offset wrt. ICRF3 position

    Parameters
    ----------
    souname: string
        source name

    Return
    ------
    dpos : a dict object
        key: "dra", "ddec"
        value: dra, ddec in mas
    """

    coordts = get_sou_ts(souname, datadir)
    dpos = calc_ts_oft(souname, coordts)
    dra = dpos["dra"]
    ddec = dpos["ddec"]
    rho = dpos["rho"]
    pa = dpos["pa"]

    coordts.add_columns([dra, ddec, rho, pa], names=["dra", "ddec", "rho", "pa"])

    return coordts


# A function deal with the source name
# --------------------------------- END --------------------------------
