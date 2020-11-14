#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_eop.py
"""
Created on Thu Mar 21 10:16:00 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Retrieve the estimates of X pole coordinate, Y pole coordinate, UT1-TAI
angle, UT1 rate, daily offsets of nutation angles as well as their formal
uncertainties and correlations from .eob file which is generated by
the program getpar.

.eob file contains series of the estimates of X pole coordinate,
Y pole coordinate, UT1-TAI angle, UT1 rate, daily offsets of nutation angles
as well as their formal uncertainties and correlations. Time tag and database
name is attached to each line. .EOB format is an extension of the IERS EOP
format.

   File contains lines of three types:
1) Comment. The first character is #. Header comments contain some information
   about solution.

2) Header. The first two symbols are blank. Header lines contain titles of the
   columns

3) Estimates.


  1    1-1    A1     ---     Usage flag
  2    3-14   F12.6  days    Modified Julian date of the TDT time tag for
                             pole coordinates and UT1
  3   16-25   A10    ---     Database name
  4   27-32   A6     ---     IVS session code (if available)
  5   34-41   F8.6   arcsec  The estimate of X pole coordinate
  6   43-50   F8.6   arcsec  The estimate of Y pole coordinate
  7   52-62   F11.7  sec     The UT1-TAI function
  8   64-71   F8.3   mas     Adjustment of the nutation in longitude angle with
                                        respect to IAU 1980 nutation expansion
  9   73-80   F8.3   mas     Adjustment of the nutation in obliquity angle with
                                        respect to IAU 1980 theory
 10   82-90   F9.6   asc/day The estimate of X pole rate
 11   92-100  F9.6   asc/day The estimate of Y pole rate
 12  102-108  F7.4   ms/day  The estimate of UT1 rate
 13  110-117  F8.6   arcsec  Formal uncertainty of X pole coordinate
 14  119-126  F8.6   arcsec  Formal uncertainty of Y pole coordinate
 15  128-136  F9.7   sec     Formal uncertainty of UT1-UTC function
 16  138-144  F7.3   mas     Formal uncertainty of nutation in longitude angle
 17  146-152  F7.3   mas     Formal uncertainty of nutation in obliquity angle
 18  154-162  F9.6   asc/day Formal uncertainty of X pole rate
 19  164-172  F9.6   asc/day Formal uncertainty of Y pole rate
 20  174-180  F7.4   asc/day Formal uncertainty of UT1 rate
 21  182-187  F6.4   --      Correlation between the estimates of X-pole
                                          positions and Y-pole position
 22  189-194  F6.4   --      Correlation between the estimates of X-pole
                                         positions and UT1-TAI angle
 23  196-201  F6.4   --      Correlation between the estimates of Y-pole
                                         positions and UT1-TAI angle
 24  203-208  F6.4   --      Correlation between the estimates of nutation in
                                         longitude and nutation in obliquity
 25  210-215  F6.4   --      Correlation between the estimates of X-pole
                                          positions and UT1 rate
 26  217-222  F6.4   --      Correlation between the estimates of Y-pole
                                         positions and UT1-TAI date
 27  224-229  F6.4   --      Correlation between the estimates of
                                         UT1-TAI angle UT1 rate
 28  231-235  F5.2   hours   Session duration
 29  237-243  F7.2   psec    Weighted root mean square of postfit residuals
 30  245-250  I6     --      Number of used observations in the session
 31  252-263  F12.6  days    Modified Julian date for nutation at TDT time
                             scale
 32  265-328  A64    --      The network configuration line. Consists of
                             two characters IVS station codes listed
                             in alphabetic order for stations that participated
                             in the experiment and supplied the data that have
                             been used in processing this experiment.

If the specific parameter was not estimated in this experiment, the field
for its value and formal uncertainty is replaced by filler: $$$$$$. The filler
takes entire field.
"""


import numpy as np
import sys
import os

from astropy.table import Table, Column
from astropy import units as u
from astropy.units import cds

__all__ = {"read_eop", "read_eob", "read_eops"}


# -----------------------------  FUNCTIONS -----------------------------
def read_eop(eop_file):
    """Retrieve from the result from .eob file.

    Parameters
    ----------
    eop_file : string
        name of .eop file

    Returns
    ----------
    t_eob : astropy.table object
        epoch_pmr : array, float
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
        lod : array, float
            UT1 rate (Length-of-Day, LOD), msec/yr
        xp_err : array, float
            formal uncertainty of X, mas
        yp_err : array, float
            formal uncertainty of Y, mas
        ut1_err : array, float
            formal uncertainty of U, msec
        dX_err : array. float
            formal uncertainty of XR, mas
        dY_err : array. float
            formal uncertainty of dY, mas
        xp_rate_err : array, float
            formal uncertainty of X-pole rate, mas/yr
        yp_rate_err : array, float
            formal uncertainty of Y-pole rate, mas/yr
        lod_err : array, float
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

    if not os.path.isfile(eop_file):
        print("Couldn't find the file", eop_file)
        sys.exit()

    t_eop = Table.read(eop_file, format="ascii.fixed_width_no_header",
                       names=["db_name", "epoch", "num_obs",
                              "xp", "xp_err",
                              "yp", "yp_err",
                              "ut1_tai", "ut1_err",
                              "xp_rate", "xp_rate_err",
                              "yp_rate", "yp_rate_err",
                              "lod", "lod_err"],
                       col_starts=[10, 33, 57, 68, 83, 98, 113, 128, 143,
                                   158, 173, 188, 203, 218, 233],
                       col_ends=[20, 49, 63, 79, 93, 109, 123, 139, 153,
                                 169, 183, 199, 213, 229, 243])

    # Add the unit information
    # 1) Time tag
    t_eop["epoch"].unit = cds.MJD

    # 2) polar motion (wobble) and rate
    t_eop["xp"].unit = u.mas
    t_eop["yp"].unit = u.mas
    t_eop["xp_err"].unit = u.uas
    t_eop["yp_err"].unit = u.uas
    t_eop["xp_rate"].unit = u.mas / u.day
    t_eop["yp_rate"].unit = u.mas / u.day
    t_eop["xp_rate_err"].unit = u.mas / u.day
    t_eop["yp_rate_err"].unit = u.mas / u.day

    # 3) UT1-TAI/UT1-UTC and rate
    t_eop["ut1_tai"].unit = u.second / 1000
    t_eop["ut1_err"].unit = u.second / 1000
    t_eop["lod"].unit = u.second / 1000
    t_eop["lod_err"].unit = u.second / 1000

    # Remove the data points which are not estimated in the solution
    mask = ((t_eop["xp_err"] != 0) & (t_eop["yp_err"] != 0)
            & (t_eop["ut1_err"] != 0))
    t_eop = Table(t_eop[mask], masked=False)

    # Sort the eob series chronologically
    t_eop.sort("epoch_pmr")

    return t_eop


def read_eob(eob_file):
    """Retrieve from the result from .eob file.

    Parameters
    ----------
    eob_file : string
        name of .eob file

    Returns
    ----------
    t_eob : astropy.table object
        epoch_pmr : array, float
            time lag, modified Julian date
        xp : array, float
            X-pole coordinate, as
        yp : array, float
            Y-pole coordinate, as
        ut1_tai : array, float
            UT1 - TAI, sec
        dX : array, float
            X coordinate of nutation offset, mas
        dY : array, float, mas
            Y coordinate of Nutation offset, mas
        xp_rate : array, float
            X-pole rate, as/yr
        yp_rate : array, float
            Y-pole rate, as/yr
        lod : array, float
            UT1 rate, msec/day
        xp_err : array, float
            formal uncertainty of X, mas
        yp_err : array, float
            formal uncertainty of Y, mas
        ut1_err : array, float
            formal uncertainty of U, msec
        dX_err : array. float
            formal uncertainty of XR, mas
        dY_err : array. float
            formal uncertainty of dY, mas
        xp_rate_err : array, float
            formal uncertainty of X-pole rate, as/yr
        yp_rate_err : array, float
            formal uncertainty of Y-pole rate, as/yr
        lod_err : array, float
            formal uncertainty of UT1 rate, msec/day
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

    if not os.path.isfile(eob_file):
        print("Couldn't find the file", eob_file)
        sys.exit()

    t_eob = Table.read(eob_file, format="ascii",
                       names=["epoch_pmr", "db_name",
                              "xp", "yp", "ut1_tai", "dX", "dY",
                              "xp_rate", "yp_rate", "lod",
                              "xp_err", "yp_err", "ut1_err",
                              "dX_err", "dY_err",
                              "xp_rate_err", "yp_rate_err", "lod_err",
                              "xp_yp_corr", "xp_ut1_corr",
                              "yp_ut1_corr", "dX_dY_corr",
                              "ut1_ut1rate_corr", "xp_ut1rate_corr",
                              "yp_ut1rate_corr",
                              "sess_dur", "wrms", "num_obs",
                              "epoch_nut", "network"])

    # Add the unit information
    # 1) Time tag
    t_eob["epoch_pmr"].unit = cds.MJD
    t_eob["epoch_nut"].unit = cds.MJD

    # 2) polar motion (wobble) and rate
    t_eob["xp"].unit = u.arcsec
    t_eob["yp"].unit = u.arcsec
    t_eob["xp_err"].unit = u.arcsec
    t_eob["yp_err"].unit = u.arcsec
    t_eob["xp_rate"].unit = u.arcsec / u.day
    t_eob["yp_rate"].unit = u.arcsec / u.day
    t_eob["xp_rate_err"].unit = u.arcsec / u.day
    t_eob["yp_rate_err"].unit = u.arcsec / u.day

    # 3) UT1-TAI/UT1-UTC and rate
    t_eob["ut1_tai"].unit = u.second
    t_eob["ut1_err"].unit = u.second
    t_eob["lod"].unit = u.second / 1000
    # t_eob["lod_err"] = t_eob["lod_err"] / 15  # arcsec -> second
    # I suppose this is s mistake in the help document
    t_eob["lod_err"].unit = u.second / 1000

    # 4) Nutation offset
    t_eob["dX"].unit = u.mas
    t_eob["dY"].unit = u.mas
    t_eob["dX_err"].unit = u.mas
    t_eob["dY_err"].unit = u.mas

    # Other information
    t_eob["sess_dur"].unit = u.hour
    t_eob["wrms"].unit = u.ps

    # Convert the unit
    t_eob["xp_err"].convert_unit_to(u.mas)
    t_eob["yp_err"] = t_eob["yp_err"].to(u.mas)
    t_eob["xp_rate_err"].convert_unit_to(u.mas / u.day)
    t_eob["yp_rate_err"].convert_unit_to(u.mas / u.day)
    t_eob["ut1_err"].convert_unit_to(u.second / 1e3)
    # t_eob["lod_err"].convert_unit_to(u.second / 1e3)

    # Remove the data points which are not estimated in the solution
    mask = ((t_eob["xp_err"] != 0) & (t_eob["yp_err"] != 0)
            & (t_eob["ut1_err"] != 0) & (t_eob["dX_err"] != 0)
            & (t_eob["dY_err"] != 0))
    t_eob = Table(t_eob[mask], masked=False)

    # Sort the eob series chronologically
    t_eob.sort("epoch_pmr")

    return t_eob


def read_eops(eops_file):
    """Retrieve from the result from .eops file.

    Parameters
    ----------
    eops_file : string
        name of .eops file

    Returns
    ----------
    t_eob : astropy.table object
        epoch_pmr : array, float
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
        lod : array, float
            UT1 rate, msec/yr
        xp_err : array, float
            formal uncertainty of X, mas
        yp_err : array, float
            formal uncertainty of Y, mas
        ut1_err : array, float
            formal uncertainty of U, msec
        dX_err : array. float
            formal uncertainty of XR, mas
        dY_err : array. float
            formal uncertainty of dY, mas
        xp_rate_err : array, float
            formal uncertainty of X-pole rate, mas/yr
        yp_rate_err : array, float
            formal uncertainty of Y-pole rate, mas/yr
        lod_err : array, float
            formal uncertainty of LOD rate, msec/yr
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

    if not os.path.isfile(eops_file):
        print("Couldn't find the file", eops_file)
        sys.exit()

    t_eops = Table.read(eops_file, format="ascii",
                        names=["epoch",
                               "xp", "yp", "ut1_tai", "dX", "dY",
                               "xp_err", "yp_err", "ut1_err",
                               "dX_err", "dY_err",
                               "db_name",
                               "xp_yp_corr", "xp_ut1_corr",
                               "yp_ut1_corr", "dX_dY_corr",
                               "num_obs", "sess_id",
                               "xp_rate", "yp_rate", "lod",
                               "nonused1", "nonused2",
                               "xp_rate_err", "yp_rate_err", "lod_err",
                               "wrms", "network"],
                        exlcude_name=["nonused1", "nonused2",
                                      "nonused3", "nonused4"])

    # Add the unit information
    # 1) Time tag
    t_eops["epoch"].unit = cds.MJD

    # 2) polar motion (wobble) and rate
    t_eops["xp"].unit = u.arcsec
    t_eops["yp"].unit = u.arcsec
    t_eops["xp_err"].unit = u.arcsec
    t_eops["yp_err"].unit = u.arcsec
    t_eops["xp_rate"].unit = u.arcsec / u.day
    t_eops["yp_rate"].unit = u.arcsec / u.day
    t_eops["xp_rate_err"].unit = u.arcsec / u.day
    t_eops["yp_rate_err"].unit = u.arcsec / u.day

    # 3) UT1-TAI/UT1-UTC and LOD
    t_eops["ut1_tai"].unit = u.second
    t_eops["ut1_err"].unit = u.second
    t_eops["lod"].unit = u.second / 1000
    # t_eops["lod_err"] = t_eops["lod_err"] / 15
    t_eops["lod_err"].unit = u.second / 1000

    # 4) Nutation offset
    t_eops["dX"].unit = u.mas
    t_eops["dY"].unit = u.mas
    t_eops["dX_err"].unit = u.mas
    t_eops["dY_err"].unit = u.mas

    # Other information
    t_eops["wrms"].unit = u.ps

    # Remove the data points which are not estimated in the solution
    mask = ((t_eops["xp_err"] != 0) & (t_eops["yp_err"] != 0)
            & (t_eops["ut1_err"] != 0) & (t_eops["dX_err"] != 0)
            & (t_eops["dY_err"] != 0))
    t_eops = Table(t_eops[mask], masked=False)

    # Sort the eob series chronologically
    t_eops.sort("epoch_pmr")

    return t_eops


# -------------------------------- MAIN --------------------------------
if __name__ == '__main__':
    read_eob(sys.argv[1])
# --------------------------------- END --------------------------------
