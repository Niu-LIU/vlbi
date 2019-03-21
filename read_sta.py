#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:52:21 2017

@author: Neo

Retrieve the estimates of positions of global stations and the formal
uncertainties of these estimates from .sta file which is generated by
the program getpar.


   .sta file contains estimates of positions of global stations and the formal
uncertainties of these estimates. The list of station positions is sorted in
alphabetic order of station names. Stations before and after episodic motions
are treated as different stations. Correlations between station positions and
velocities are also written.

   File contains lines of four types:

1) Comment. The first character is #. Header comment contain the full name of
   the spool file.

2) Cartesian components of the vector of station position. The first
   8 characters of this line are STA_GCX:

   Field   Format Units Meaning
   1-8     A8     --    record type identifier: STA_GCX:
   11-25   A15    --    station name. Station name consist of 8-letters station
                        acronym and 6-letter epoch in format yymmdd. Epoch
                        is attached to the name only if the station had
                        episodic motion. Fields between the last letter of the
                        station name and the first letter of epoch are filled
                        by _. If the station didn't have episodic name then
                        the space after the last letter of the station name is
                        left blank.
   28-29   A2     --    component identifier. One of "X:", "Y:" or "Z:"
   31-45   F15.2  mm    value of X-component of station position.
   50-59   F10.3  mm    formal uncertainty of X-component of station position.
   65-79   F15.2  mm    value of Y-component of station position.
   84-93   F10.3  mm    formal uncertainty of Y-component of station position.
   99-113  F15.2  mm    value of Z-component of station position.
   118-127 F10.3  mm    formal uncertainty of Z-component of station position.
   139-145 I7     --    the number of observations of this station used in
                        solution.
   156-162 I7     --    total number of observations of this station.
   174-178 I5     --    the number of sessions with this station used in
                        solution.
   189-193 I5     --    total number of sessions of this station.
   205-214 A10    --    the date of the first session with this station used
                        in solution. format: yyyy.mm.dd (as integer numbers)
   226-235 A10    --    the date of the last session with this station used
                        in solution. format: yyyy.mm.dd (as integer numbers)

3) Local topocentric components of the vector of station position: Up, East,
   North. The first 8 characters of this line are STA_GCU:

   Field   Format Units Meaning
   1-8     A8     --    record type identifier: STA_GCU:
   11-25   A15    --    station name. Station name consist of 8-letters station
                        acronym and 6-letter epoch in format yymmdd. Epoch
                        is attached to the name only if the station had
                        episodic motion. Fields between the last letter of the
                        station name and the first letter of epoch are filled
                        by _. If the station didn't have episodic name then
                        the space after the last letter of the station name is
                        left blank.
   28-29   A2     --    component identifier. One of "U:", "E:" or "N:"
   31-45   F15.2  mm    value of U-component of station position.
   50-59   F10.3  mm    formal uncertainty of U-component of station position.
   65-79   F15.2  mm    value of E-component of station position.
   84-93   F10.3  mm    formal uncertainty of E-component of station position.
   99-113  F15.2  mm    value of N-component of station position.
   118-127 F10.3  mm    formal uncertainty of N-component of station position.

4) Correlations between station positions and velocities. Correlation matrix
   is defined as the matrix of 6x6 in the upper triangle representation without
   the main diagonal which. Elements in the columns or rows of the matrix are
   in the order: X-position, Y-position, Z-position, X-velocity, Y-velocity,
   Z-velocity.

   1-8     A8     --    record type identifier: STA_CRL:
   11-25   A15    --    station name. Station name consist of 8-letters station
                        acronym and 6-letter epoch in format yymmdd. Epoch
                        is attached to the name only if the station had
                        episodic motion. Fields between the last letter of the
                        station name and the first letter of epoch are filled
                        by _. If the station didn't have episodic name then
                        the space after the last letter of the station name is
                        left blank.
   31-36   F6.3   d/l   Correlation between X-position and Y-position
   38-43   F6.3   d/l   Correlation between X-position and Z-position
   45-50   F6.3   d/l   Correlation between Y-position and Z-position
   52-57   F6.3   d/l   Correlation between X-position and X-velocity
   59-64   F6.3   d/l   Correlation between Y-position and X-velocity
   66-71   F6.3   d/l   Correlation between Z-position and X-velocity
   73-78   F6.3   d/l   Correlation between X-position and Y-velocity
   80-85   F6.3   d/l   Correlation between Y-position and Y-velocity
   87-92   F6.3   d/l   Correlation between Z-position and Y-velocity
   94-99   F6.3   d/l   Correlation between X-velocity and Y-velocity
   101-106 F6.3   d/l   Correlation between X-position and Z-velocity
   108-113 F6.3   d/l   Correlation between Y-position and Z-velocity
   115-120 F6.3   d/l   Correlation between Z-position and Z-velocity
   122-127 F6.3   d/l   Correlation between X-velocity and Z-velocity
   129-134 F6.3   d/l   Correlation between Y-velocity and Z-velocity

"""

from astropy.table import Table, join
from astropy import units as u
import numpy as np
import os
import sys
from convert_func import date2jyear


__all__ = ["read_sta"]


# ------------------------------  FUNCTIONS  ---------------------------
def read_sta_gcx(gcx_file):
    '''Retrieve cartesian components of the vector of station positions.

    Parameters
    ----------
    gcx_file : string
        name of data file

    Returns
    ----------
    sta_gcx : astropy.table object
     |__
        station : string
            name of station
        xp : array, float
            X component of position
        yp : array, float
            Y component of position
        zp : array, float
            Z component of position
        xp_err : array, float
            formal uncertainty of X component
        yp_err : array, float
            formal uncertainty of Y component
        zp_err : array, float
            formal uncertainty of Z component
        used_obs : array,
            Number of used observations of this source
        total_obs : array, int
            Total number of observations of this source
        used_sess : array, int
            Number of used sessions for this source
        yotal_sess : array, int
            Total number of sessions for this source
        beg_date : array, float
            Epoch of the first observation
        end_date : array, float
            Epoch of the last observation
    '''

    # Check if the input file exists
    if not os.path.exists(gcx_file):
        print("Couldn't find the input file", gcx_file)
        sys.exit()

    sta_gcx = Table.read(gcx_file,
                         format="ascii.fixed_width_no_header",
                         names=["station",
                                "xp", "xp_err", "yp", "yp_err", "zp", "zp_err",
                                "used_obs", "total_obs", "used_sess",
                                "total_sess", "beg_date", "end_date"],
                         col_starts=[10, 30, 49, 64, 83, 98, 117, 138, 155,
                                     173, 188, 204, 225],
                         col_ends=[25, 45, 59, 79, 93, 113, 127, 145, 162,
                                   178, 193, 214, 235])

    # Filling the missing values
    sta_gcx["beg_date"] = sta_gcx["beg_date"].filled(" "*10)
    sta_gcx["end_date"] = sta_gcx["end_date"].filled(" "*10)

    # # Date transformation
    beg_date_list = [date2jyear(beg_datei)
                     for beg_datei in sta_gcx["beg_date"]]
    sta_gcx["beg_date"] = beg_date_list

    end_date_list = [date2jyear(end_datei)
                     for end_datei in sta_gcx["end_date"]]
    sta_gcx["end_date"] = end_date_list

    # Add information for units
    sta_gcx["xp"].unit = u.m / 1000
    sta_gcx["yp"].unit = u.m / 1000
    sta_gcx["zp"].unit = u.m / 1000
    sta_gcx["xp_err"].unit = u.m / 1000
    sta_gcx["yp_err"].unit = u.m / 1000
    sta_gcx["zp_err"].unit = u.m / 1000
    sta_gcx["beg_date"].unit = u.yr
    sta_gcx["end_date"].unit = u.yr

    return sta_gcx


def read_sta_gcu(gcu_file):
    '''Retrievel local topocentric components of station positions.

    Parameters
    ----------
    gcu_file : string
        name of data file

    Returns
    ----------
    sta_gcx : astropy.table object
     |__
        station : string
            name of station
        up : array, float
            U component
        ep : array, float
            E component
        np : array, float
            N component
        up_err : array, float
            formal uncertainty of U component
        ep_err : array, float
            formal uncertainty of E component
        np_err : array, float
            formal uncertainty of N component
    '''

    # Check if the input file exists
    if not os.path.exists(gcu_file):
        print("Couldn't find the input file", gcu_file)
        sys.exit()

    sta_gcu = Table.read(gcu_file,
                         format="ascii.fixed_width_no_header",
                         names=["station", "up", "up_err", "ep", "ep_err",
                                "np", "np_err"],
                         col_starts=[10, 30, 49, 64, 83, 98, 117],
                         col_ends=[25, 45, 59, 79, 93, 113, 127])

    # Add information for units
    sta_gcu["up"].unit = u.m / 1000
    sta_gcu["ep"].unit = u.m / 1000
    sta_gcu["np"].unit = u.m / 1000
    sta_gcu["up_err"].unit = u.m / 1000
    sta_gcu["ep_err"].unit = u.m / 1000
    sta_gcu["np_err"].unit = u.m / 1000

    return sta_gcu


def read_sta_crl(crl_file):
    '''Retrievel local topocentric components of station positions.

    Parameters
    ----------
    crl_file : string
        name of data file

    Returns
    ----------
    sta_gcx : astropy.table object
     |__
        station : string
            name of station
        xp_yp_corr : array, float
            Correlation between X-position and Y-position
        xp_zp_corr : array, float
            Correlation between X-position and Z-position
        yp_zp_corr : array, float
            Correlation between Y-position and Z-position
        xp_xv_corr : array, float
            Correlation between X-position and X-velocity
        yp_xv_corr : array, float
            Correlation between Y-position and X-velocity
        zp_xv_corr : array, float
            Correlation between Z-position and X-velocity
        xp_yv_corr : array, float
            Correlation between X-position and Y-velocity
        yp_yv_corr : array, float
            Correlation between Y-position and Y-velocity
        zp_yv_corr : array, float
            Correlation between Z-position and Y-velocity
        xv_yv_corr : array, float
            Correlation between X-velocity and Y-velocity
        xp_zv_corr : array, float
            Correlation between X-position and Z-velocity
        yp_zv_corr : array, float
            Correlation between Y-position and Z-velocity
        zp_zv_corr : array, float
            Correlation between Z-position and Z-velocity
        xv_zv_corr : array, float
            Correlation between X-velocity and Z-velocity
        yv_zv_corr : array, float
            Correlation between Y-velocity and Z-velocity
    '''

    # Check if the input file exists
    if not os.path.exists(crl_file):
        print("Couldn't find the input file", crl_file)
        sys.exit()

    sta_crl = Table.read(crl_file,
                         format="ascii.fixed_width_no_header",
                         names=["station",
                                "xp_yp_corr", "xp_zp_corr", "yp_zp_corr",
                                "xp_xv_corr", "yp_xv_corr", "zp_xv_corr",
                                "xp_yv_corr", "yp_yv_corr", "zp_yv_corr",
                                "xv_yv_corr", "xp_zv_corr", "yp_zv_corr",
                                "zp_zv_corr", "xv_zv_corr", "yv_zv_corr"],
                         col_starts=[10, 30, 37, 44, 51, 58, 65, 73, 79, 86,
                                     93, 100, 107, 114, 121, 128],
                         col_ends=[25, 36, 43, 50, 57, 64, 71, 78, 85, 92,
                                   99, 106, 113, 120, 127, 134])

    return sta_crl


def parse_sta_file(in_file, out_file="temp.out", keyword=None):
    '''printer the lines containing keyword in .sta file

    Parameters
    ----------
    in_file : string
        .sta file
    out_file : string
        output file, default is None
    keyword : string
        keyword for grep. It could be set as "STA_CRL", "STA_GCU", or
        "STA_GCX".
    '''

    keyword_list = ["STA_GCX", "STA_GCU", "STA_CRL"]

    if keyword is None or not keyword in keyword_list:
        print("Valid keyword!")
        sys.exit()

    # Check if the input file exists
    if not os.path.exists(in_file):
        print("Couldn't find the input file", in_file)
        sys.exit()

    # Assume that the operation system is an Unix- or Linux-like where
    # the commond grep is supported.
    print("Generate a temporary file", out_file, "of type", keyword)
    os.system("grep %s %s > %s" % (keyword, in_file, out_file))
    print("Done!")

    # Read the data
    print("Read", out_file)
    if keyword is "STA_GCX":
        data_table = read_sta_gcx(out_file)
    elif keyword is "STA_GCU":
        data_table = read_sta_gcu(out_file)
    else:
        data_table = read_sta_crl(out_file)
    print("Done!")

    # Delete the temporary file
    print("Delete", out_file)
    os.system("rm %s" % out_file)
    print("Done!")

    return data_table


def read_sta(sta_file):
    '''Retrieve the result from .sta file.

    Parameters
    ----------
    sta_file : string
        name of .sta file

    Returns
    ----------
    sta_pos : astropy.table object
     |__
        station : string
            name of station
        xp : array, float
            X component
        yp : array, float
            Y component
        zp : array, float
            Z component
        xp_err : array, float
            formal uncertainty of X component
        yp_err : array, float
            formal uncertainty of Y component
        zp_err : array, float
            formal uncertainty of Z component
        up : array, float
            U component
        ep : array, float
            E component
        np : array, float
            N component
        up_err : array, float
            formal uncertainty of U component
        ep_err : array, float
            formal uncertainty of E component
        np_err : array, float
            formal uncertainty of N component
        xp_yp_corr : array, float
            Correlation between X-position and Y-position
        xp_zp_corr : array, float
            Correlation between X-position and Z-position
        yp_zp_corr : array, float
            Correlation between Y-position and Z-position
        xp_xv_corr : array, float
            Correlation between X-position and X-velocity
        yp_xv_corr : array, float
            Correlation between Y-position and X-velocity
        zp_xv_corr : array, float
            Correlation between Z-position and X-velocity
        xp_yv_corr : array, float
            Correlation between X-position and Y-velocity
        yp_yv_corr : array, float
            Correlation between Y-position and Y-velocity
        zp_yv_corr : array, float
            Correlation between Z-position and Y-velocity
        xv_yv_corr : array, float
            Correlation between X-velocity and Y-velocity
        xp_zv_corr : array, float
            Correlation between X-position and Z-velocity
        yp_zv_corr : array, float
            Correlation between Y-position and Z-velocity
        zp_zv_corr : array, float
            Correlation between Z-position and Z-velocity
        xv_zv_corr : array, float
            Correlation between X-velocity and Z-velocity
        yv_zv_corr : array, float
            Correlation between Y-velocity and Z-velocity
        used_obs : array,
            Number of used observations of this source
        total_obs : array, int
            Total number of observations of this source
        used_sess : array, int
            Number of used sessions for this source
        yotal_sess : array, int
            Total number of sessions for this source
        beg_date : array, float
            Epoch of the first observation
        end_date : array, float
            Epoch of the last observation
    '''

    # Original data contains 3 types of data, which will be read repsectively
    # 1) STA_GCX (Cartesian components of the vector of station position)
    sta_gcx = parse_sta_file(sta_file, keyword="STA_GCX")

    # 2) STA_GCU (Local topocentric components of station position)
    sta_gcu = parse_sta_file(sta_file, keyword="STA_GCU")

    # 3) STA_CRL (Correlations between station positions and velocities)
    sta_crl = parse_sta_file(sta_file, keyword="STA_CRL")

    # Merge three tables into one major table
    table0 = join(sta_gcx, sta_gcu, keys="station", join_type="outer")
    sta_pos = join(table0, sta_crl, keys="station", join_type="outer")

    return sta_pos


if __name__ == '__main__':
    print("No thing to do!")
    pass
# ------------------------------ END -----------------------------------
