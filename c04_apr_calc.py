#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: c04_apr_calc.py
"""
Created on Sun Dec 10 18:37:34 2017

@author: Neo(liuniu@smail.nju.edu.cn)

Calculate the apriori EOP based C04 series.

History
14/06/2018 : print differences of EOP and Nutation parameters into a
             single file
"""

from astropy.table import Table, Column
from astropy import units as u
from astropy.units import cds
import numpy as np
from scipy.interpolate import CubicSpline
import sys
# My module
from .read_c04 import read_c04_array


__all__ = {"read_c04", "calc_c04_apr"}


# -----------------------------  FUNCTIONS -----------------------------
def read_c04(c04_file="/Users/Neo/tmp/git/Niu-LIU/vlbi/aux_files/"
             "eopc04_IAU2000.62-now.txt"):
    """Read EOP from C04 series.

    Parameters
    ----------
    c04_file : string
        c04 data file

    Returns
    ----------
    t_c04 : astropy.table object
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


def read_c04(c04_file="/Users/Neo/tmp/git/Niu-LIU/vlbi/aux_files/"
             "eopc04_IAU2000.62-now.txt"):
    """Read EOP from C04 series.

    Parameters
    ----------
    c04_file : string
        c04 data file

    Returns
    ----------
    t_c04 : astropy.table object
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


def read_c04_array(c04_file="/Users/Neo/tmp/git/Niu-LIU/vlbi/aux_files/"
                   "eopc04_IAU2000.62-now.txt"):
    """Fetch the C04 series.

    Parameters
    ----------
    c04_file : string
                path and name of C04 data file.

    Returns
    ----------
    epoch_mjd : array, float
        epoch in modified Julian date
    xp : array, float
        xp position of CIP in ITRS, as
    yp : array, float
        yp position of CIP in ITRS, as
    ut : array, float
        UT1 - UTC, ms
    # LOD : array, float
    #     length of day, ms
    dX : array, float
        xp component of CPO, mas
    dY : array, float
        yp component of CPO, mas
    """

    epoch_mjd, xp, yp, ut = np.genfromtxt(c04_file, skip_header=14,
                                          usecols=np.arange(3, 7), unpack=True)
    dX, dY = np.genfromtxt(c04_file, skip_header=14,
                           usecols=(8, 9), unpack=True)

    # Convert unit from arc-sec to mas to be consistent with .eob file
    dX, dY = dX * 1000, dY * 1000

    return epoch_mjd, xp, yp, ut, dX, dY


def read_c04_apr(c04_apr_file):
    """Read EOP from C04 series.

    Parameters
    ----------
    c04_apr_file : string
        c04 a priori file

    Returns
    ----------
    t_c04 : astropy.table object
    """

    if not os.path.isfile(c04_apr_file):
        print("Couldn't find the file", c04_apr_file)
        sys.exit()

    t_c04 = Table.read(c04_apr_file, format="ascii",
                       names=["epoch",
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


def cubic_spline(x, xs, ys):
    '''For series (xi,yi), get y(x) using cubic spline.


    Parameters
    ----------
    xs : array, float
        time series of X component
    ys : array, float
        time series of Y component
    x : float
        X position need to be interpolated

    Returns
    ---------
    Interpolated Y value.
    '''

    cs = CubicSpline(xs, ys)

    return cs(x)


def interpolate_pmr(epoch_pmr, epoch_c04, xp_c04, yp_c04, ut_c04):
    """Get interpolated EOP at a certain epoch.


    Parameters
    ----------
    epoch_pmre : float
        epoch to be interpolated
    poch_c04 : array, float
        20-point epoch series
    xp_sub : array, float
        20-point xp series
    yp_sub : array, float
        20-point yp series
    ut_sub : array, float
        20-point UT1-UTC series

    Returns
    ----------
    xp_apr : float
        interpolated xp value
    yp_apr : float
        interpolated yp value
    dut_apr : float
        interpolated UT1-UTC value
    """

    xp_apr = np.zeros_like(epoch_pmr)
    yp_apr = np.zeros_like(epoch_pmr)
    ut_apr = np.zeros_like(epoch_pmr)

    epo_ind = np.searchsorted(epoch_c04, epoch_pmr)

    for i, ind in enumerate(epo_ind):
        if ind == 0 or ind >= epoch_c04.size:
            # normally it won't happen!!
            print("The epoch %f was too early or too late"
                  " for the C04 series." % epoch_pmr[i])
            sys.exit()

        elif ind < 9 or epoch_c04.size - ind < 10:
            # In this case we will use less datapoints.
            pnum = np.min(ind, epoch_c04.size - ind)

        else:
            pnum = 10

        epoch_c04_sub = epoch_c04[ind - pnum + 1: ind + pnum + 1]

        # EOP sub-series
        xp_sub = xp_c04[ind - pnum + 1: ind + pnum + 1]
        yp_sub = yp_c04[ind - pnum + 1: ind + pnum + 1]
        ut_sub = ut_c04[ind - pnum + 1: ind + pnum + 1]

        # Cubic interpolation
        xp_apr[i] = cubic_spline(epoch_pmr[i], epoch_c04_sub, xp_sub)
        yp_apr[i] = cubic_spline(epoch_pmr[i], epoch_c04_sub, xp_sub)
        ut_apr[i] = cubic_spline(epoch_pmr[i], epoch_c04_sub, xp_sub)

    return xp_apr, yp_apr, dut_apr


def interpolate_nut(epoch_pmr, epoch_c04, dX_c04, dY_c04):
    '''Get interpolated EOP at a certain epoch.


    Parameters
    ----------
    epo : float
        epoch to be interpolated
    mjd20 : array, float
        20-point epoch series, Julian day
    dX20 : array, float
        20-point dX series, mas
    dY20 : array, float
        20-point dY series, mas

    Returns
    ----------
    dX_apr : float
        interpolated dX value, mas
    dY_apr : float
        interpolated dY value, mas
    '''

    # dX
    dX_apr = cubic_spline(epo, mjd20, dX20)
    # dY
    dY_apr = cubic_spline(epo, mjd20, dY20)

    return dX_apr, dY_apr


def calc_c04_apr(t_eob, ofile=None):
    """Use Cubic Spline Interpolation to calculate the EOP.


    For the target eopch, We use the 10 successive points in front of
    and behind it in the time series (total 20 points) to do the Cubic
    Spline Interpolation.

    Parameters
    ----------
    t_eob : astropy.table object
        EOP estimates from .eob file
    """

    print("# ---------- BEGIN ----------")

    epoch_pmr = np.array(t_eob["epoch_pmr"])
    epoch_nut = np.array(t_eob["epoch_nut"])

    # Read C04 data
    [epoch_c04, xp_c04, yp_c04, ut_c04, dX_c04, dY_c04] = read_c04_array()

    # Output
    if ofile is None:
        ofile = "temp.eob_c04_apr"
    print("# Output file: %s" % (ofile))

    fout = open(ofile, "w")
    # header
    print("# EOP and Nutation a priori values based C04:\n"
          "# Epoch  xp   yp   UT1-UTC  ddX  ddY\n"
          "# mjd    as   as   second   mas  mas",
          file=fout)

    # Interpolation
    print("\n# Begin to interpolate!")

    # For polar motion and UT1
    # epo_ind = np.searchsorted(epoch_c04, epoch_pmr)

    # xp_apr = []
    # yp_apr = []
    # dut_apr = []

    # for epoch_pmri, ind in zip(epoch_pmr, epo_ind):
    #     if ind == 0 or ind >= epoch_c04.size:
    #         # normally it won't happen!!
    #         print("The epoch %f was too early or too late"
    #               " for the C04 series." % epochi)
    #         sys.exit()

    #     elif ind < 9 or epoch_c04.size - ind < 10:
    #         # In this case we will use less datapoints.
    #         pnum = np.min(ind, epoch_c04.size - ind)

    #     else:
    #         pnum = 10

    #     mjd20 = epoch_c04[ind - pnum + 1: ind + pnum + 1]

    #     # EOP
    #     xp_20 = xp[ind - pnum + 1: ind + pnum + 1]
    #     yp_20 = yp[ind - pnum + 1: ind + pnum + 1]
    #     ut_20 = ut[ind - pnum + 1: ind + pnum + 1]

    #     xp_apri, yp_apri, dut_apri = interpolate_pmr(
    #         epochi, mjd20, xp_20, yp_20, ut_20)

    # Nutation
    dX20 = dX[ind - pnum + 1: ind + pnum + 1]
    dY20 = dY[ind - pnum + 1: ind + pnum + 1]

    dX_apr, dY_apr = interpolate_nut(
        epochi, mjd20, dX20, dY20)

    print("%f  %+8.3f  %+8.3f  %+8.3f  %+8.3f  %+8.3f" %
          (epochi, xp_apr, yp_apr, dut_apr, dX_apr, dY_apr), file=fout)

    fout.close()

    print("# ----------  END  ----------")


# --------------------------------- MAIN -------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 2:
        calc_c04_apr(sys.argv[1])
    else:
        print("Input error!")
# --------------------------------- END --------------------------------
