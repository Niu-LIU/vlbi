#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_sou.py
"""
Created on Sat Sep 15 14:38:46 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, unique, Column
import astropy.units as u
from astropy.units import cds
from astropy.coordinates import Angle
import numpy as np
import os
import sys
import time

# My modules
from my_progs.vlbi.convert_func import RA_conv, DC_conv, date2mjd
from my_progs.catalog.pos_err import pos_err_calc

__all__ = ["read_sou", "write_cat", "rewrite_sou", "read_cat"]


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
    """

    if not os.path.isfile(sou_file):
        sys.exit()

    # print("  Remove these sources with 0 observation used in Solve.")
    # os.system("./sou_elim %s" % sou_file)

    t_sou0 = Table.read(sou_file, format="ascii.fixed_width_no_header",
                        names=["ivs_name", "ra", "ra_err",
                               "dec", "dec_err", "ra_dec_corr",
                               "obs_used", "obs_tot", "ses_used",
                               "ses_tot", "date_beg", "date_end"],
                        col_starts=(10, 24, 45, 61, 82, 98, 117,
                                    132, 150, 164, 181, 202),
                        col_ends=(18, 41, 55, 78, 92, 104, 122,
                                  139, 155, 170, 191, 212))

    mask = (t_sou0["obs_used"] != 0)
    t_sou1 = Table(t_sou0[mask], masked=False)

    # Table of source name (IERS/IVS name)
    t_name = Table.read("/Users/Neo/Tools/packages/calc_solve/mk5-opa/"
                        "save_files/source.names",
                        format="ascii.fixed_width_no_header",
                        names=["ivs_name", "iers_name"],
                        col_starts=(0, 22),
                        col_ends=(8, 30))

    t_name = unique(t_name, keys="ivs_name")

    # Join these two tables
    t_sou = join(t_sou1, t_name, keys="ivs_name", join_type="left")

    # Change the order of columns
    # Copy the columns
    iers_name = t_sou["iers_name"]

    # Remove the existing columns
    t_sou.remove_column("iers_name")

    # insert these columns
    t_sou.add_column(iers_name, 1)

    # Fill the masked(missing) value
    t_sou["iers_name"].fill_value = "-" * 8
    t_sou["iers_name"] = t_sou["iers_name"].filled()

    # convert string into float for RA, Decl. and observing epoch
    Nsize = len(t_sou)
    ra = np.empty(shape=(Nsize,), dtype=float)
    dec = np.empty(shape=(Nsize,), dtype=float)
    date_beg = np.empty(shape=(Nsize,), dtype=float)
    date_end = np.empty(shape=(Nsize,), dtype=float)

    for i, t_soui in enumerate(t_sou):
        ra[i] = RA_conv(t_soui["ra"])
        dec[i] = DC_conv(t_soui["dec"])
        date_beg[i] = date2mjd(t_soui["date_beg"])
        date_end[i] = date2mjd(t_soui["date_end"])

    # replace original columns with new columns
    t_sou["ra"] = ra
    t_sou["dec"] = dec
    t_sou["date_beg"] = date_beg
    t_sou["date_end"] = date_end

    # unit
    t_sou["ra"].unit = u.deg
    t_sou["dec"].unit = u.deg
    t_sou["dec_err"].unit = u.mas
    t_sou["date_beg"].unit = cds.MJD
    t_sou["date_end"].unit = cds.MJD

    # Multipliy the formal error in R.A. by a factor of cos(Decl.)
    factor = np.cos(Angle(t_sou["dec"]).radian)
    t_sou["ra_err"] = t_sou["ra_err"] * factor
    t_sou["ra_err"].unit = u.mas

    # Calculate the positional error
    pos_err_sx = pos_err_calc(t_sou["ra_err"], t_sou["dec_err"],
                              t_sou["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    t_sou.add_column(pos_err_sx, name="pos_err", index=6)
    t_sou["pos_err"].unit = u.mas

    return t_sou


def write_cat(t_sou, sou_file, cat_file):
    """
    Parameters
    ----------
    t_sou : astropy.table object
    cat_file : string
        the full path of .sou file

    """

    # New columns
    date_mean = Column((t_sou["date_beg"] + t_sou["date_end"]) / 2.,
                       name="date_mean")

    # insert these columns
    t_sou.add_column(date_mean, 10)

    # Add comments
    t_sou.meta["comments"] = [
        "VLBI Celestial Reference Frame Solution quasarpm-180510a",
        " Columns  Units   Meaning",
        "    1     --      IVS designation",
        "    2     --      IERS designation",
        "    3     deg     Right ascension",
        "    4     mas     Formal uncertainty of the right ascension "
        "(*cos(Dec))",
        "    5     deg     Declination",
        "    6     mas     Formal uncertainty of the declination",
        "    7     --      Correlation between right ascension and "
        "declination",
        "    8     --      Number of delays",
        "    9     --      Number of sessions",
        "   10     days    Average epoch of observation (MJD)",
        "   11     days    First epoch of observation (MJD)",
        "   12     days    Last epoch of observation (MJD)",
        " Created date: %s." % time.strftime("%d/%m/%Y", time.localtime())]

    t_sou.write(cat_file, format="ascii.fixed_width_no_header",
                delimiter="  ",
                exclude_names=["obs_tot", "ses_tot"],
                formats={"ivs_name": "%8s", "iers_name": "%8s",
                         "ra": "%14.10f", "dec": "%+14.10f",
                         "ra_err": "%10.4f", "dec_err": "%10.4f",
                         "ra_dec_corr": "%+7.4f"},
                overwrite=True)
    print("Rewrite results in %s into %s : Done!" % (sou_file, cat_file))


def rewrite_sou(sou_file):
    """Rewrite the source position into .cat file

    Parameters
    ----------
    sou_file : string
        the full path of .sou file

    Return
    ------
    t_sou : astropy.table object
    """

    # Read .sou file to get source positions
    t_sou = read_sou(sou_file)

    # write the source position into a .cat file
    cat_file = "%s.cat" % sou_file[:-4]

    write_cat(t_sou, sou_file, cat_file)


def read_cat(cat_file):
    """Read radio source positions

    Parameters
    ----------
    cat_file : string
        the full path of .cat file

    Return
    ------
    t_sou : astropy.table object
    """

    if not os.path.isfile(cat_file):
        sys.exit()

    t_sou0 = Table.read(cat_file, format="ascii",
                        names=["ivs_name", "iers_name", "ra",
                               "dec", "ra_err", "dec_err", "ra_dec_corr",
                               "mean_epo", "beg_epo", "end_epo",
                               "num_sess", "num_del", "num_delrate", "flag"])

    mask = (t_sou0["num_del"] != 0)
    t_sou = Table(t_sou0[mask], masked=False)

    # convert string into float for RA, Decl. and observing epoch
    Nsize = len(t_sou)

    # unit
    t_sou["ra"].unit = u.deg
    t_sou["dec"].unit = u.deg
    t_sou["ra_err"].unit = u.mas
    t_sou["dec_err"].unit = u.mas
    t_sou["mean_epo"].unit = cds.MJD
    t_sou["beg_epo"].unit = cds.MJD
    t_sou["end_epo"].unit = cds.MJD

    # Calculate the positional error
    pos_err_sx = pos_err_calc(t_sou["ra_err"], t_sou["dec_err"],
                              t_sou["ra_dec_corr"])

    # Add the semi-major axis of error ellipse to the table
    t_sou.add_column(pos_err_sx, name="pos_err", index=7)

    return t_sou


# -------------------------------- MAIN --------------------------------
if __name__ == '__main__':
    rewrite_sou(sys.argv[1])
# --------------------------------- END --------------------------------
"""
   Field   Format Units     Meaning
   1-8     A8     --        record type identifier: EOP_LOC:
   11-20   A10    --        database name with leading dollar sign
   23-25   I3     --        database version number
   34-49   A16    calend    EOP time tag in Solve format: YYYY.MM.DD-hh:mm
                            Time scale is not defined. Adjustments are at TDB
                            time scale, a priori EOP are at unknown time scale.
   58-63   I6     --        number of observation used for getting these EOP
                            estimates.
   69-79   F11.4  mas       estimate of X-pole coordinate
   84-93   F10.2  muas      formal uncertainty of X-pole coordinate
   99-109  F11.4  mas       estimate of Y-pole coordinate
   114-123 F10.2  muas      formal uncertainty of Y-pole coordinate
   129-139 F11.4  msec      estimates of UT1-TAI
   144-153 F10.2  musec     formal uncertainty of UT1-TAI
   159-169 F11.4  mas/day   estimates of X pole rate
   174-183 F10.2  muas/day  formal uncertainties of X pole rate
   189-199 F11.4  msec/day  estimates of Y pole rate
   204-213 F10.2  msec/day  formal uncertainties of Y pole rate
   219-229 F11.4  msec/day  estimates of UT1-TAI rate
   234-243 F10.2  musec/day formal uncertainties of UT1-TAI rate
   249-259 F11.4  ms/day**2 estimates of UT1-TAI acceleration
   264-273 F10.2  ms/day**2 formal uncertainties of UT1-TAI acceleration
"""
