#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_sou.py
"""
Created on Sat Sep 15 14:38:46 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from astropy.table import Table, join, unique, Column
import astropy.units as u
from astropy.units import cds
from astropy.coordinates import Angle
import numpy as np
import os
import sys
import time

# My modules
# from .convert_func import RA_conv, DC_conv, date2mjd
# from my_progs.catalog.pos_err import pos_err_calc
from convert_func import RA_conv, DC_conv, date2mjd
from pos_err import pos_err_calc

__all__ = ["read_sou", "read_cat"]


# -----------------------------  FUNCTIONS -----------------------------
def read_sou(sou_file):
    """Read radio source positions

    Parameters
    ----------
    sou_file : string
        the full path of .sou file

    Return
    ------
    t_sou : astropy.table object
        |
        -- ivs_name : str
            IVS source name
        -- ra : float
            right ascension (degree)
        -- ra_err
            formal uncertainty in RA (mas)
        -- dec
            declination (degree)
        -- dec_err
            formal uncertainty in Dec. (mas)
        -- ra_dec_corr
            correlation coefficient between RA and Dec.
        -- pos_err
            ellipse semi-major axis of positional error (mas)
        -- used_obs
            number of observations used in the solution
        -- total_obs
            number of total observations of this source
        -- used_sess
            number of sessions used in the solution
        -- total_sess
            number of total sessions of this source
        -- beg_date
            epoch of first observation (MJD)
        -- end_date
            epoch of last observation (MJD)
    """

    if not os.path.isfile(sou_file):
        sys.exit()

    t_sou = Table.read(sou_file, format="ascii.fixed_width_no_header",
                       names=("ivs_name", "ra", "ra_err",
                              "dec", "dec_err", "ra_dec_corr",
                              "used_obs", "total_obs", "used_sess",
                              "total_sess", "beg_date", "end_date"),
                       col_starts=(10, 24, 45, 61, 82, 98, 117,
                                   132, 150, 164, 181, 202),
                       col_ends=(18, 41, 55, 78, 92, 104, 122,
                                 139, 155, 170, 191, 212))

    ra_dec_table = Table.read(sou_file, format="ascii.fixed_width_no_header",
                              names=["ra", "dec", "used_obs", ],
                              col_starts=[20, 57, 117], col_ends=[41, 78, 122])

    # Remove source with 0 observation used in th solution
    mask = (t_sou["used_obs"] != 0)
    t_sou = Table(t_sou[mask], masked=False)
    mask = (ra_dec_table["used_obs"] != 0)
    ra_dec_table = Table(ra_dec_table[mask], masked=False)

    # convert string into float for RA, Decl. and observing epoch
    Nsize = len(t_sou)
    ra = np.empty(shape=(Nsize,), dtype=float)
    dec = np.empty(shape=(Nsize,), dtype=float)
    date_beg = np.empty(shape=(Nsize,), dtype=float)
    date_end = np.empty(shape=(Nsize,), dtype=float)

    for i, ra_dec_tablei in enumerate(ra_dec_table):
        ra[i] = RA_conv(ra_dec_tablei["ra"])
        dec[i] = DC_conv(ra_dec_tablei["dec"])

    for i, t_soui in enumerate(t_sou):
        date_beg[i] = date2mjd(t_soui["beg_date"])
        date_end[i] = date2mjd(t_soui["end_date"])

    # replace original columns with new columns
    t_sou["ra"] = ra
    t_sou["dec"] = dec
    t_sou["beg_date"] = date_beg
    t_sou["end_date"] = date_end

    # Add unit information
    t_sou["ra"].unit = u.deg
    t_sou["dec"].unit = u.deg
    t_sou["ra_err"].unit = u.mas
    t_sou["dec_err"].unit = u.mas
    t_sou["beg_date"].unit = cds.MJD
    t_sou["end_date"].unit = cds.MJD

    # Multipliy the formal error in R.A. by a factor of cos(Decl.)
    factor = np.cos(Angle(t_sou["dec"]).radian)
    t_sou["ra_err"] = t_sou["ra_err"] * factor

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(t_sou["ra_err"], t_sou["dec_err"],
                           t_sou["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    t_sou.add_column(pos_err, name="pos_err", index=7)
    t_sou["pos_err"].unit = u.mas

    return t_sou


def read_cat(cat_file):
    """Read radio source positions

    Parameters
    ----------
    cat_file : string
        the full path of .cat file

    Return
    ------
    t_sou : astropy.table object
        |
        -- ivs_name : str
            IVS source name
        -- ra : float
            right ascension (degree)
        -- ra_err
            formal uncertainty in RA (mas)
        -- dec
            declination (degree)
        -- dec_err
            formal uncertainty in Dec. (mas)
        -- ra_dec_corr
            correlation coefficient between RA and Dec.
        -- pos_err
            ellipse semi-major axis of positional error (mas)
        -- used_obs
            number of observations used in the solution
        -- total_obs
            number of total observations of this source
        -- used_sess
            number of sessions used in the solution
        -- total_sess
            number of total sessions of this source
        -- beg_date
            epoch of first observation (MJD)
        -- end_date
            epoch of last observation (MJD)
    """

    if not os.path.isfile(cat_file):
        sys.exit()

    t_sou = Table.read(cat_file, format="ascii",
                       names=["ivs_name", "iers_name", "ra",
                              "dec", "ra_err", "dec_err", "ra_dec_corr",
                              "mean_epo", "beg_epo", "end_epo",
                              "num_sess", "num_del", "num_delrate", "flag"])

    mask = (t_sou["num_del"] != 0)
    t_sou = Table(t_sou[mask], masked=False)

    # unit
    t_sou["ra"].unit = u.deg
    t_sou["dec"].unit = u.deg
    t_sou["ra_err"].unit = u.mas
    t_sou["dec_err"].unit = u.mas
    t_sou["mean_epo"].unit = cds.MJD
    t_sou["beg_epo"].unit = cds.MJD
    t_sou["end_epo"].unit = cds.MJD

    # Calculate the semi-major axis of error ellipse
    pos_err = pos_err_calc(t_sou["ra_err"], t_sou["dec_err"],
                           t_sou["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    t_sou.add_column(pos_err, name="pos_err", index=7)
    t_sou["pos_err"].unit = u.mas

    return t_sou


# -------------------------------- MAIN --------------------------------
if __name__ == '__main__':
    rewrite_sou(sys.argv[1])
# --------------------------------- END --------------------------------
