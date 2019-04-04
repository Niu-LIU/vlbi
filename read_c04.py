#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_c04.py
"""
Created on Wed Apr  3 15:06:53 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, Column
from astropy import units as u
from astropy.units import cds
import numpy as np
import sys
import os


# -----------------------------  FUNCTIONS -----------------------------
def read_c04(c04_file="/Users/Neo/tmp/git/Niu-LIU/vlbi/aux_files/"
             "eopc04_IAU2000.62-now.txt"):
    """Read EOP from C04 series.

    Parameters
    ----------
    c04_file : string
        name of .eop file

    Returns
    ----------
    t_eob : astropy.table object
        epoch : array, float
            time lag, modified Julian date
        xp : array, float
            X-pole coordinate, mas
        yp : array, float
            Y-pole coordinate, mas
        ut1_tai : array, float
            UT1 - TAI, msec
        dX : array, float
            X coordinate of nutation offset, mas
        dY : array, float, mas
            Y coordinate of Nutation offset, mas
        xp_rate : array, float
            X-pole rate, mas/yr
        yp_rate : array, float
            Y-pole rate, mas/yr
        ut1_rate : array, float
            UT1 rate, msec/yr
        xp_err : array, float
            formal uncertainty of X, mas
        yp_err : array, float
            formal uncertainty of Y, mas
        dut1_err : array, float
            formal uncertainty of U, msec
        dX_err : array. float
            formal uncertainty of XR, mas
        dY_err : array. float
            formal uncertainty of dY, mas
        xp_rate_err : array, float
            formal uncertainty of X-pole rate, mas/yr
        yp_rate_err : array, float
            formal uncertainty of Y-pole rate, mas/yr
        dut_rate_err : array, float
            formal uncertainty of UT1 rate, msec/yr
        xp_yp_corr : array, float
            correlation between xp and yp
        xp_ut1_corr : array, float
            correlation between xp and UT1
        yp_ut1_corr : array, float
            correlation between yp and UT1
        dX_dY_corr : array, float
            correlation between dX and dY
        ut1_ut1rate_corr : array, float
            correlation between UT1 and UT1 rate
        xp_ut1rate_err : array, float
            correlation between xp and UT1 rate
        yp_ut1rate_corr : array, float
            correlation between yp and UT1 rate
        sess_dur : array, float
            Session duration, hours
        wrms : array, float
            wrms of postfit residuals, ps
        num_obs : array, float
            number of used observations in the session
        epoch_nut : array, float
            mjd for nutation at TDT time scale
        network : array, string
            network configuration line
    """

    if not os.path.isfile(c04_file):
        print("Couldn't find the file", c04_file)
        sys.exit()

    t_c04 = Table.read(c04_file, format="ascii",
                       names=["year", "month", "date", "epoch",
                              "xp", "yp", "ut1_utc", "lod", "dX", "dY",
                              "xp_err", "yp_err", "dut1_err",
                              "lod_err", "dX_err", "dY_err"])

    # Add the unit information
    # 1) Time tag
    t_c04["epoch"].unit = cds.MJD

    # 2) polar motion (wobble) and rate
    t_c04["xp"].unit = u.arcsec
    t_c04["yp"].unit = u.arcsec
    t_c04["xp_err"].unit = u.arcsec
    t_c04["yp_err"].unit = u.arcsec

    # 3) UT1-UTC and LOD
    t_c04["ut1_utc"].unit = u.second
    t_c04["dut1_err"].unit = u.second
    t_c04["lod"].unit = u.second
    t_c04["lod_err"].unit = u.second

    # 4) polar motion (wobble) and rate
    t_c04["dX"].unit = u.arcsec
    t_c04["dY"].unit = u.arcsec
    t_c04["dX_err"].unit = u.arcsec
    t_c04["dY_err"].unit = u.arcsec

    # Remove the data points which are not estimated in the solution
    mask = ((t_c04["xp_err"] != 0) & (t_c04["yp_err"] != 0)
            & (t_c04["dut1_err"] != 0) & (t_c04["lod_err"] != 0)
            & (t_c04["dX_err"] != 0) & (t_c04["dY_err"] != 0))

    t_c04 = Table(t_c04[mask], masked=False)

    return t_c04


def read_c04_array(C04_file="aux_files/eopc04_IAU2000.62-now.txt"):
    '''Fetch the C04 series.

    Parameters
    ----------
    C04_file : string
                path and name of C04 data file.

    Returns
    ----------
    mjd : array, float
        epoch in modified Julian date
    Xp : array, float
        Xp position of CIP in ITRS, mas
    Yp : array, float
        Yp position of CIP in ITRS, mas
    U : array, float
        UT1 - UTC, ms
    # LOD : array, float
    #     length of day, ms
    dX : array, float
        Xp component of CPO, mas
    dY : array, float
        Yp component of CPO, mas
    XpErr : array, float
        formal uncertainty of Xp, mas
    YpErr : array, float
        formal uncertainty of Yp, mas
    UErr : array, float
        formal uncertainty of U, ms
    # LODErr : array, float
    #     formal uncertainty of LOD, ms
    dXErr : array, float
        formal uncertainty of dX, mas
    dYErr : array, float
        formal uncertainty of dY, mas
    '''

    mjd, Xp, Yp, U = np.genfromtxt(C04_file, skip_header=14,
                                   usecols=np.arange(3, 7), unpack=True)
    dX, dY = np.genfromtxt(C04_file, skip_header=14,
                           usecols=(8, 9), unpack=True)
    XpErr, YpErr, UErr = np.genfromtxt(C04_file, skip_header=14,
                                       usecols=np.arange(10, 13),
                                       unpack=True)
    dXErr, dYErr = np.genfromtxt(C04_file, skip_header=14,
                                 usecols=(14, 15),
                                 unpack=True)

    # arc-sec -->  mas, second --> ms
    Xp, Yp, U = Xp * 1000, Yp * 1000, U * 1000
    XpErr, YpErr, UErr = XpErr * 1000, YpErr * 1000, UErr * 1000

    dX, dY = dX * 1000, dY * 1000
    dXErr, dYErr = dXErr * 1000, dYErr * 1000

    return mjd, Xp, Yp, U, dX, dY, XpErr, YpErr, UErr, dXErr, dYErr


# -------------------------------- MAIN --------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 2:
        read_c04(sys.argv[1])
    else:
        print("Input error!")
# --------------------------------- END --------------------------------
