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

from astropy.table import Table, Column, join
from astropy import units as u
from astropy.units import cds
import numpy as np
import os
from scipy.interpolate import CubicSpline
import sys
# My module
from .delta_tai_utc import delta_tai_utc_calc


__all__ = {"read_c04", "calc_c04_apr", "calc_c04_offset"}


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
    epoch_pmr : float
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

    insert_ind = np.searchsorted(epoch_c04, epoch_pmr)

    for i, ind in enumerate(insert_ind):

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
        yp_apr[i] = cubic_spline(epoch_pmr[i], epoch_c04_sub, yp_sub)
        ut_apr[i] = cubic_spline(epoch_pmr[i], epoch_c04_sub, ut_sub)

        # In C04, ut_apr is UT1-UTC while it is UT1-TAI in .eob file,
        # so we need to calculate IAI - UTC and then convert ut_apr to UT1-UTC
        delta_tai_utc = delta_tai_utc_calc(epoch_pmr[i])

        ut_apr[i] -= delta_tai_utc

    return xp_apr, yp_apr, ut_apr


def interpolate_nut(epoch_nut, epoch_c04, dX_c04, dY_c04):
    """Get interpolated Nutation offset at a certain epoch.


    Parameters
    ----------
    epoch_nut : float
        epoch to be interpolated
    poch_c04 : array, float
        20-point epoch series
    dX_sub : array, float
        20-point dX series
    dY_sub : array, float
        20-point dY series

    Returns
    ----------
    dX_apr : float
        interpolated xp value
    dY_apr : float
        interpolated yp value
    """

    dX_apr = np.zeros_like(epoch_nut)
    dY_apr = np.zeros_like(epoch_nut)

    insert_ind = np.searchsorted(epoch_c04, epoch_nut)

    for i, ind in enumerate(insert_ind):
        if ind == 0 or ind >= epoch_c04.size:
            # normally it won't happen!!
            print("The epoch %f was too early or too late"
                  " for the C04 series." % epoch_nut[i])
            sys.exit()

        elif ind < 9 or epoch_c04.size - ind < 10:
            # In this case we will use less datapoints.
            pnum = np.min(ind, epoch_c04.size - ind)

        else:
            pnum = 10

        epoch_c04_sub = epoch_c04[ind - pnum + 1: ind + pnum + 1]

        # EOP sub-series
        dX_sub = dX_c04[ind - pnum + 1: ind + pnum + 1]
        dY_sub = dY_c04[ind - pnum + 1: ind + pnum + 1]

        # Cubic interpolation
        dX_apr[i] = cubic_spline(epoch_nut[i], epoch_c04_sub, dX_sub)
        dY_apr[i] = cubic_spline(epoch_nut[i], epoch_c04_sub, dY_sub)

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

    Returns
    -------
    t_c04_apr : astropy.table object
        C04 a priori values for EOPs
    """

    epoch_pmr = np.array(t_eob["epoch_pmr"])
    epoch_nut = np.array(t_eob["epoch_nut"])

    # Read C04 data
    [epoch_c04, xp_c04, yp_c04, ut_c04, dX_c04, dY_c04] = read_c04_array()

    # Interpolation

    # For polar motion and UT1
    xp_apr, yp_apr, dut_apr = interpolate_pmr(
        epoch_pmr, epoch_c04, xp_c04, yp_c04, ut_c04)

    # For nutation offset
    dX_apr, dY_apr = interpolate_nut(epoch_nut, epoch_c04, dX_c04, dY_c04)

    # Creat a table to store the a priori information
    t_c04_apr = Table(
        [t_eob["db_name"], xp_apr, yp_apr, dut_apr, dX_apr, dY_apr],
        names=["db_name", "xp_c04_apr", "yp_c04_apr", "ut1_tai_c04_apr",
               "dX_c04_apr", "dY_c04_apr"])

    # Add the unit information
    # 1) polar motion (wobble) and rate
    t_c04_apr["xp_c04_apr"].unit = u.arcsec
    t_c04_apr["yp_c04_apr"].unit = u.arcsec

    # 2) UT1-TAI/UT1-UTC and rate
    t_c04_apr["ut1_tai_c04_apr"].unit = u.second

    # 3) Nutation offset
    t_c04_apr["dX_c04_apr"].unit = u.mas
    t_c04_apr["dY_c04_apr"].unit = u.mas

    # Output
    if ofile is not None:
        print("# Output file: %s" % (ofile))

        fout = open(ofile, "w")
        # header
        print("# EOP and Nutation a priori values based C04:\n"
              "# Epoch  xp   yp   UT1-UTC  ddX  ddY\n"
              "# mjd    as   as   second   mas  mas",
              file=fout)

        for (epochi, xp_apri, yp_apri, dut_apri, dX_apri, dY_apri) in zip(
                epoch_pmr, xp_apr, yp_apr, dut_apr, dX_apr, dY_apr):
            print("%f  %+8.3f  %+8.3f  %+8.3f  %+8.3f  %+8.3f" %
                  (epochi, xp_apr, yp_apr, dut_apr, dX_apr, dY_apr), file=fout)

        fout.close()

    return t_c04_apr


def calc_c04_offset(t_eob, ofile=None):
    """Calculate the EOP offset wrt. the C04 series.

    Parameters
    ----------
    t_eob : astropy.table object
        EOP estimates from .eob file

    Returns
    -------
    t_eob : astropy.table object
        EOP estimates from .eob file
    """

    # C04 a priori values
    t_c04_apr = calc_c04_apr(t_eob)

    # Cross-match two tables by sessoion identifier
    t_eob_apr = join(t_eob, t_c04_apr, keys="db_name")

    dxp = t_eob_apr["xp_c04_apr"] - t_eob_apr["xp"]
    dyp = t_eob_apr["yp_c04_apr"] - t_eob_apr["yp"]
    dut = t_eob_apr["ut1_tai_c04_apr"] - t_eob_apr["ut1_tai"]
    ddX = t_eob_apr["dX_c04_apr"] - t_eob_apr["dX"]
    ddY = t_eob_apr["dY_c04_apr"] - t_eob_apr["dY"]

    # Polar motion
    dxp.convert_unit_to(u.uas)
    dyp.convert_unit_to(u.uas)

    # UT1-UTC (s -> uas)
    dut = dut * 15e6
    dut.unit = u.uas

    # Nutation offset
    ddX.convert_unit_to(u.uas)
    ddY.convert_unit_to(u.uas)

    t_eob_apr.add_columns([dxp, dyp, dut, ddX, ddY],
                          names=["dxp_c04", "dyp_c04", "dut_c04",
                                 "ddX_c04", "ddY_c04"])

    # Output
    if ofile is not None:
        print("# Output file: %s" % (ofile))

        fout = open(ofile, "w")
        # header
        print("# EOP and Nutation difference wrt C04\n"
              "# Epoch  dXp  err  dYp  err  dUT  err  ddX  err  ddY  err\n"
              "# mjd    mas  mas  mas  mas  mas  mas  mas  mas  mas  mas",
              file=fout)

        # EOP offset
        for (i, mjd) in enumerate(t_eob_apr["epoch_pmr"]):
            print("%13.6f" % mjd, "  %+10.3f  %7.3f" * 5 %
                  (dxp[i], t_eob_apr["xp_err"][i],
                   dyp[i], t_eob_apr["yp_err"][i],
                   dut[i], t_eob_apr["dut1_err"][i],
                   ddX[i], t_eob_apr["dX_err"][i],
                   ddY[i], t_eob_apr["dY_err"][i]), file=fout)

        fout.close()

    print("# ----------  END  ----------")

    return t_eob_apr


# --------------------------------- MAIN -------------------------------
if __name__ == "__main__":
    print("Nothing to do")
# --------------------------------- END --------------------------------
