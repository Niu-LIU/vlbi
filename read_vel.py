#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:55:55 2017

@author: Neo

Retrieve the estimates of velocities of global stations and the formal
uncertainties of these estimates from .vel file which is generated by
the program getpar.

   .vel file contains values estimates of velocities of global stations and the
formal uncertainties of these estimates. The list of the estimates is sorted in
alphabetic order of station names. Stations before and after episodic motions
are treated as different stations. Correlations between station positions and
velocities are also written.

   File contains lines of three types:

1) Comment. The first character is #. Header comment contain the full name of
   the spool file.

2) Cartesian components of the vector of station velocity. The first
   8 characters of this line are STA_GVX:

   Field   Format Units Meaning
   1-8     A8     --    record type identifier: STA_GVX:
   11-18   A8     --    station name.
   24-32   F9.2   mm/yr value of X-component of station velocity.
   37-44   F8.3   mm/yr formal uncertainty of X-component of station velocity.
   50-58   F9.2   mm/yr value of Y-component of station velocity.
   63-70   F8.3   mm/yr formal uncertainty of Y-component of station velocity.
   76-84   F9.2   mm/yr value of Z-component of station velocity.
   89-96   F8.3   mm/yr formal uncertainty of Z-component of station velocity.

3) Local topocentric components of the vector of station velocity: Up, East,
   North. The first 8 characters of this line are STA_GVU:

   Field   Format Units Meaning
   1-8     A8     --    record type identifier: STA_GVU:
   11-18   A8     --    station name.
   24-32   F9.2   mm/yr value of U-component of station velocity.
   37-44   F8.3   mm/yr formal uncertainty of U-component of station velocity.
   50-58   F9.2   mm/yr value of E-component of station velocity.
   63-70   F8.3   mm/yr formal uncertainty of E-component of station velocity.
   76-84   F9.2   mm/yr value of N-component of station velocity.
   89-96   F8.3   mm/yr formal uncertainty of N-component of station velocity.
"""


from astropy.table import Table, join
from astropy import units as u
import numpy as np
import os
import sys
# My module
from .read_sta import read_sta


__all__ = ["read_vel"]


# ------------------------------  FUNCTIONS  ---------------------------
def read_sta_gvx(gvx_file):
    '''Retrieve cartesian components of station velocities.

    Parameters
    ----------
    gvx_file : string
        name of data file

    Returns
    ----------
    sta_gvx : astropy.table object
     |__
        station : string
            name of station
        xv : array, float
            X component of station velocity
        xv_err : array, float
            formal uncertainty of X component of station velocity
        yv : array, float
            Y component of station velocity
        yv_err : array, float
            formal uncertainty of Y component of station velocity
        zv : array, float
            Z component of station velocity
        zv_err : array, float
            formal uncertainty of Z component of station velocity
    '''

    # Check if the input file exists
    if not os.path.exists(gvx_file):
        print("Couldn't find the input file", gvx_file)
        sys.exit()

    sta_gvx = Table.read(gvx_file,
                         format="ascii.fixed_width_no_header",
                         names=["station", "xv", "xv_err", "yv", "yv_err",
                                "zv", "zv_err"],
                         col_starts=[10, 23, 36, 49, 62, 75, 88],
                         col_ends=[18, 32, 44, 58, 70, 84, 96])

    # Add information for units
    sta_gvx["xv"].unit = u.m / 1000 / u.yr
    sta_gvx["yv"].unit = u.m / 1000 / u.yr
    sta_gvx["zv"].unit = u.m / 1000 / u.yr
    sta_gvx["xv_err"].unit = u.m / 1000 / u.yr
    sta_gvx["yv_err"].unit = u.m / 1000 / u.yr
    sta_gvx["zv_err"].unit = u.m / 1000 / u.yr

    return sta_gvx


def read_sta_gvu(gvu_file):
    '''Retrieve local topocentric components of station velocities.

    Parameters
    ----------
    gvu_file : string
        name of data file

    Returns
    ----------
    sta_gvu : astropy.table object
     |__
        station : string
            name of station
        uv : array, float
            U component of station velocity
        uv_err : array, float
            formal uncertainty of U component of station velocity
        ev : array, float
            E component of station velocity
        ev_err : array, float
            formal uncertainty of E component of station velocity
        nv : array, float
            N component of station velocity
        nv_err : array, float
            formal uncertainty of N component of station velocity
    '''

    # Check if the input file exists
    if not os.path.exists(gvu_file):
        print("Couldn't find the input file", gvu_file)
        sys.exit()

    sta_gvu = Table.read(gvu_file,
                         format="ascii.fixed_width_no_header",
                         names=["station", "uv", "uv_err", "ev", "ev_err",
                                "nv", "nv_err"],
                         col_starts=[10, 23, 36, 49, 62, 75, 88],
                         col_ends=[18, 32, 44, 58, 70, 84, 96])

    # Add information for units
    sta_gvu["uv"].unit = u.m / 1000 / u.yr
    sta_gvu["ev"].unit = u.m / 1000 / u.yr
    sta_gvu["nv"].unit = u.m / 1000 / u.yr
    sta_gvu["uv_err"].unit = u.m / 1000 / u.yr
    sta_gvu["ev_err"].unit = u.m / 1000 / u.yr
    sta_gvu["nv_err"].unit = u.m / 1000 / u.yr

    return sta_gvu


def parse_vel_file(in_file, out_file="temp.out", keyword=None):
    '''printer the lines containing keyword in .vel file

    Parameters
    ----------
    in_file : string
        .sta file
    out_file : string
        output file, default is None
    keyword : string
        keyword for grep. It could be set as "STA_GVU" or "STA_GVX".
    '''

    keyword_list = ["STA_GVX", "STA_GVU"]

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
    if keyword is "STA_GVX":
        data_table = read_sta_gvx(out_file)
    else:
        data_table = read_sta_gvu(out_file)
    print("Done!")

    # Delete the temporary file
    print("Delete", out_file)
    os.system("rm %s" % out_file)
    print("Done!")

    return data_table


def read_vel(vel_file, get_corr=False):
    '''Retrieve the result from .lso file.

    Parameters
    ----------
    vel_file : string
        name of .vel file

    Returns
    ----------
    sta_vel : astropy.table object
     |__
        station : string
            name of station
        xv : array, float
            X component of station velocity
        yv : array, float
            Y component of station velocity
        zv : array, float
            Z component of station velocity
        xv_err : array, float
            formal uncertainty of X component of station velocity
        yv_err : array, float
            formal uncertainty of Y component of station velocity
        zv_err : array, float
            formal uncertainty of Z component of station velocity
        uv : array, float
            U component of station velocity
        ev : array, float
            E component of station velocity
        nv : array, float
            N component of station velocity
        uv_err : array, float
            formal uncertainty of U component of station velocity
        ev_err : array, float
            formal uncertainty of E component of station velocity
        nv_err : array, float
            formal uncertainty of N component of station velocity
    '''

    # Original data contains 2 types of data, which will be read repsectively
    # 1) STA_GVX (Cartesian components of the vector of station velocities)
    sta_gvx = parse_vel_file(vel_file, keyword="STA_GVX")

    # 2) STA_GVU (Local topocentric components of station velocities)
    sta_gvu = parse_vel_file(vel_file, keyword="STA_GVU")

    # Merge two tables into one major table
    sta_vel = join(sta_gvx, sta_gvu, keys="station", join_type="outer")

    # Do not forget the correlation between estimates.
    if get_corr:
        sta_corr = read_sta("%s.sta" % vel_file[:-4])
        sta_corr.keep_columns(
            ["station", "xv_yv_corr",  "xv_zv_corr", "yv_zv_corr"])

        sta_vel = join(sta_gvx, sta_corr, keys="station", join_type="left")

    return sta_vel


# ----------------------------- MAIN -----------------------------------
if __name__ == '__main__':
    print("Nothing to do!")
# ------------------------------ END -----------------------------------
